"""Codon alignment with optimized gap placement.

Provides a pure-Python/Cython codon alignment implementation and a CLI
command that supports switching between the Python backend and the
Rust-accelerated backend (``codon_alignment_rust``).

Algorithm overview
------------------
1. **Gap gathering** — nearby gaps within ``min_gap_distance`` are merged
   into single windows and moved to the centre of each window.
2. **Codon grouping** — the reference and sequence are split into codon
   triples, and consecutive gap-containing codons are grouped.
3. **Gap placement optimisation** — for each gap group the algorithm
   searches every codon-aligned position (step 3) using a centre-expand
   order starting from the original gap location.  At each candidate
   position a combined IUPAC + BLOSUM62 score is computed.  Precomputed
   other-side amino acids avoid redundant translation work.
4. **Finalisation** — gaps inside sequence codons are moved to the codon
   end so downstream consumers always see gaps trailing the bases.

Optimisations over the original implementation (D1-D4):
- **D1** precomputed other-side amino acids,
- **D2** fast-path when no gaps exist,
- **D4** centre-expand search order with leftmost-index tie-breaking.
"""

import re
from collections.abc import Iterable
from itertools import chain, groupby
from typing import Any

import click
import cython  # type: ignore

from ..cli import cli
from ..models import NAPosition, RefSeqPair, Sequence
from ..processor import Processor, intermediate_processor
from ..utils import find_codon_trim_slice, group_by_codons
from ..utils.blosum62 import blosum62_score
from ..utils.codonutils import translate_codons
from ..utils.iupac import iupac_score

NOGAP: int = 0b00
REFGAP: int = 0b01
SEQGAP: int = 0b10

LEFT: int = 0b00
RIGHT: int = 0b01

GAP_PLACEMENT_SCORE_PATTERN: re.Pattern = re.compile(
    r'^(\d+)(?:/(\d+))?(ins|del):(-?\d+)$'
)

CodonPair = tuple[
    int,  # refpos0
    tuple[
        list[NAPosition],  # refcodon
        list[NAPosition],  # poscodon
    ],
]


# -----------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def extend_codons_until_gap(
    ref_codons: list[list[NAPosition]],
    seq_codons: list[list[NAPosition]],
    direction: int,
) -> tuple[list[list[NAPosition]], list[list[NAPosition]], int]:
    """Extend codon lists in *direction* until a gap is encountered.

    Args:
        ref_codons: Reference codon lists to extend from.
        seq_codons: Sequence codon lists to extend from.
        direction: ``LEFT`` or ``RIGHT``.

    Returns:
        A 3-tuple ``(ref_codons, seq_codons, count)`` where *count* is
        the number of gap-free codons that were included.
    """
    if direction == LEFT:
        ref_codons.reverse()
        seq_codons.reverse()
    idx: int
    refcd: list[NAPosition]
    seqcd: list[NAPosition]
    endidx: int = len(ref_codons)
    for idx, (refcd, seqcd) in enumerate(
        zip(ref_codons, seq_codons, strict=False)
    ):
        broken: bool = False
        for na in chain(refcd, seqcd):
            if na.is_gap:
                endidx = idx
                broken = True
                break
        if broken:
            break
    ref_codons = ref_codons[:endidx]
    seq_codons = seq_codons[:endidx]
    if direction == LEFT:
        ref_codons.reverse()
        seq_codons.reverse()
    return ref_codons, seq_codons, len(ref_codons)


@cython.cfunc
@cython.inline
@cython.returns(list)
def find_windows_with_gap(
    refnas: list[NAPosition], seqnas: list[NAPosition], min_gap_distance: int
) -> list[slice]:
    """Return slices covering contiguous gap regions.

    Gaps separated by more than *min_gap_distance* non-gap positions
    are placed in separate windows.
    """
    refna: NAPosition
    seqna: NAPosition
    first_gap_idx: int = -1
    last_gap_idx: int = -1
    windows: list[slice] = []
    for idx, (refna, seqna) in enumerate(zip(refnas, seqnas, strict=False)):
        if not refna.is_gap and not seqna.is_gap:
            continue
        if first_gap_idx == -1:
            first_gap_idx = last_gap_idx = idx
        elif idx - last_gap_idx > min_gap_distance:
            windows.append(slice(first_gap_idx, last_gap_idx + 1))
            first_gap_idx = last_gap_idx = idx
        else:  # idx - first_gap_idx <= na_window_size
            last_gap_idx = idx
    if first_gap_idx > -1:
        windows.append(slice(first_gap_idx, last_gap_idx + 1))
    return windows


@cython.cfunc
@cython.inline
def find_first_gap(nas: list[NAPosition]) -> int:
    """Return the index of the first gap in *nas*, or ``-1``."""
    idx: int
    na: NAPosition
    for idx, na in enumerate(nas):
        if na.is_gap:
            return idx
    return -1


@cython.cfunc
@cython.inline
@cython.returns(list)
def move_gap_to_codon_end(codons: list[list[NAPosition]]) -> list[list[NAPosition]]:
    """Move gaps to the end of each codon triple."""
    new_codons: list[list[NAPosition]] = []
    for codon in codons:
        new_codons.append(
            [na for na in codon if not na.is_gap] + [na for na in codon if na.is_gap]
        )
    return new_codons


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def separate_gaps_from_nas(
    nas: list[NAPosition],
) -> tuple[list[NAPosition], list[NAPosition]]:
    """Partition *nas* into ``(non_gaps, gaps)``."""
    na: NAPosition
    nongaps: list[NAPosition] = []
    gaps: list[NAPosition] = []
    for na in nas:
        if na.is_gap:
            gaps.append(na)
        else:
            nongaps.append(na)
    return nongaps, gaps


@cython.cfunc
@cython.inline
@cython.returns(list)
def remove_n_gaps(nas: list[NAPosition], n_gaps: int) -> list[NAPosition]:
    """Remove the first *n_gaps* gap positions from *nas*."""
    na: NAPosition
    result: list[NAPosition] = []
    for na in nas:
        if na.is_gap and n_gaps > 0:
            n_gaps -= 1
        else:
            result.append(na)
    return result


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def remove_redundant_gaps(
    refnas: list[NAPosition], seqnas: list[NAPosition]
) -> tuple[list[NAPosition], list[NAPosition]]:
    """Remove gaps that appear in both ref and seq simultaneously.

    When both sequences carry gaps in the same window the minimum
    overlap is stripped from each so that gaps exist on only one side.
    """
    n_gaps: int = min(NAPosition.count_gaps(refnas), NAPosition.count_gaps(seqnas))
    if n_gaps:
        refnas = remove_n_gaps(refnas, n_gaps)
        seqnas = remove_n_gaps(seqnas, n_gaps)
    return refnas, seqnas


# -----------------------------------------------------------------------
# D1 — Optimised scoring: precompute other_aas, avoid redundant work
# -----------------------------------------------------------------------


@cython.cfunc
@cython.inline
def calc_match_score_precomputed(
    mynas: list[NAPosition],
    othernas: list[NAPosition],
    other_aas: list[bytes],
    base_score: float,
) -> float:
    """Compute combined IUPAC + BLOSUM62 match score.

    Re-uses *other_aas* (precomputed amino-acid translations of the
    other side) to avoid redundant ``translate_codons`` calls.

    Args:
        mynas: Candidate gap-inserted positions.
        othernas: The other (fixed) side positions.
        other_aas: Precomputed amino-acid list for *othernas*.
        base_score: Starting score (typically ``-gaplen``).

    Returns:
        The total alignment score for this candidate placement.
    """
    myna: NAPosition
    otherna: NAPosition
    myaa: bytes
    otheraa: bytes
    myaas: list[bytes] = translate_codons(mynas)
    score: float = base_score
    for myna, otherna in zip(mynas, othernas, strict=False):
        score += iupac_score(myna.notation, otherna.notation)
    for myaa, otheraa in zip(myaas, other_aas, strict=False):
        score += blosum62_score(myaa, otheraa)
    return score


# -----------------------------------------------------------------------
# D4 — Centre-expand search order
# -----------------------------------------------------------------------


@cython.cfunc
@cython.inline
@cython.returns(list)
def center_expand_positions(
    center: int, scanstart: int, mynas_len: int, step: int
) -> list[int]:
    """Generate search positions in centre-expand order.

    Starting from *center* (snapped to the nearest *step*-aligned
    position), positions are emitted outward in alternating
    lower/higher order until all valid indices are covered.

    Args:
        center: Original gap index to expand from.
        scanstart: Minimum valid index (3 for REFGAP, 0 for SEQGAP).
        mynas_len: Length of the non-gap portion of the sequence.
        step: Codon step size (always 3).

    Returns:
        Ordered list of candidate gap-insertion indices.
    """
    positions: list[int] = []
    # Snap center to nearest valid step-aligned position
    center = max(scanstart, ((center + step // 2) // step) * step)
    if scanstart <= center <= mynas_len:
        positions.append(center)

    offset: int = step
    while True:
        added: bool = False
        lo: int = center - offset
        hi: int = center + offset
        if lo >= scanstart:
            positions.append(lo)
            added = True
        if hi <= mynas_len and hi != lo:
            positions.append(hi)
            added = True
        if not added:
            break
        offset += step

    return positions


# -----------------------------------------------------------------------
# D1 + D4 — Optimised find_best_matches
# -----------------------------------------------------------------------


@cython.cfunc
@cython.inline
@cython.returns(list)
def find_best_matches(
    mynas: list[NAPosition],
    othernas: list[NAPosition],
    bp1_indices: set[int],
    gap_type: int,
    gap_placement_score: dict[tuple[int, int], int],
    is_start: bool,
    is_end: bool,
) -> list[NAPosition]:
    """Find the optimal gap insertion position.

    Slides the gap block across every codon-aligned position using
    centre-expand order (D4) and scores each candidate with
    precomputed other-side amino acids (D1).

    Tie-breaking priority (highest first):
    1. Codon-boundary insertion (``bp1_indices``) — +1 bonus.
    2. Original gap position preserved.
    3. Leftmost index (``-idx``).

    Args:
        mynas: The side that contains the gap.
        othernas: The other (fixed) side.
        bp1_indices: Indices at reading-frame codon boundaries.
        gap_type: ``REFGAP`` or ``SEQGAP``.
        gap_placement_score: ``{(napos, gaplen): score}`` bonus map.
        is_start: Whether this window touches the sequence start.
        is_end: Whether this window touches the sequence end.

    Returns:
        The re-ordered *mynas* with the gap at the best position.
    """
    idx: int
    score: tuple[float, int, int]
    mygap: list[NAPosition]
    test_mynas: list[NAPosition]
    orig_gapidx: int = find_first_gap(mynas)
    mynas, mygap = separate_gaps_from_nas(mynas)
    gaplen: int = len(mygap)
    max_score: tuple[float, int, int] | None = None
    best_mynas: list[NAPosition] | None = None
    scanstart: int = 3 if gap_type == REFGAP else 0
    mynas_len: int = len(mynas)

    # D1: Precompute other-side amino acids (they never change)
    other_aas: list[bytes] = translate_codons(othernas)

    # D4: Search from original gap position outward
    positions: list[int] = center_expand_positions(orig_gapidx, scanstart, mynas_len, 3)

    for idx in positions:
        napos: int
        test_mynas = mynas[::]
        test_mynas[idx:idx] = mygap
        base_score: float = float(-gaplen)
        if (is_start and idx == 0) or (is_end and idx + 3 > mynas_len):
            base_score = 0.0
        # D1: Use precomputed other_aas
        score_val: float = calc_match_score_precomputed(
            test_mynas, othernas, other_aas, base_score
        )
        napos = mynas[idx - 1].pos if gap_type == REFGAP else othernas[idx].pos
        if (napos, gaplen) in gap_placement_score:
            score_val += gap_placement_score[(napos, gaplen)]
        elif (napos, 0) in gap_placement_score:
            score_val += gap_placement_score[(napos, 0)]

        if idx in bp1_indices:
            # reward gaps inserted between codons
            score = (score_val + 1, 2, -idx)
        elif idx == orig_gapidx:
            # respect the original gapidx if it's already one of the best
            score = (score_val, 1, -idx)
        else:
            score = (score_val, 0, -idx)
        if max_score is None or score > max_score:
            max_score = score
            best_mynas = test_mynas
    if best_mynas is None:
        # fallback to mynas, if no best match is found
        best_mynas = mynas
    return best_mynas


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def paired_find_best_matches(
    refnas: list[NAPosition],
    seqnas: list[NAPosition],
    gap_type: int,
    gap_placement_score: dict[int, dict[tuple[int, int], int]],
    is_seq_start: bool,
    is_seq_end: bool,
) -> tuple[list[NAPosition], list[NAPosition]]:
    """Dispatch gap optimisation to the correct side.

    Computes reading-frame boundary indices (``bp1_indices``) from the
    reference, then calls ``find_best_matches`` on whichever side
    carries the gap.

    Args:
        refnas: Reference positions (may contain gaps if REFGAP).
        seqnas: Sequence positions (may contain gaps if SEQGAP).
        gap_type: ``REFGAP`` or ``SEQGAP``.
        gap_placement_score: Full ``{gap_type: {(pos, len): score}}`` map.
        is_seq_start: True if window touches sequence start.
        is_seq_end: True if window touches sequence end.

    Returns:
        Updated ``(refnas, seqnas)`` with the gap repositioned.
    """
    idx: int
    na: NAPosition
    bp1_indices: set[int] = set()
    bp: int = 0
    for idx, na in enumerate(refnas):
        if na.is_gap:
            continue
        bp = (bp + 1) % 3
        if bp == 1:
            bp1_indices.add(idx)

    if gap_type == REFGAP:
        refnas = find_best_matches(
            refnas,
            seqnas,
            bp1_indices,
            gap_type,
            gap_placement_score[gap_type],
            # for REFGAPs, ending gaps also have penalty
            False,
            False,
        )
    elif gap_type == SEQGAP:
        seqnas = find_best_matches(
            seqnas,
            refnas,
            bp1_indices,
            gap_type,
            gap_placement_score[gap_type],
            is_seq_start,
            is_seq_end,
        )
    return refnas, seqnas


@cython.ccall
def codon_pairs_group_key(
    cdpair: tuple[int, tuple[list[NAPosition], list[NAPosition]]],
) -> int:
    """Return the gap type for a ``(index, (refcodon, seqcodon))`` pair."""
    refcd: list[NAPosition]
    seqcd: list[NAPosition]
    _, (refcd, seqcd) = cdpair
    if NAPosition.any_has_gap(refcd):
        return REFGAP
    if NAPosition.any_has_gap(seqcd):
        return SEQGAP
    return NOGAP


@cython.cfunc
@cython.inline
@cython.returns(list)
def move_gaps_to_center(nas: list[NAPosition]) -> list[NAPosition]:
    """Move all gaps in *nas* to the centre of the non-gap bases."""
    new_nas: list[NAPosition]
    gaps: list[NAPosition]
    new_nas, gaps = separate_gaps_from_nas(nas)
    center_idx: int = len(new_nas) // 2
    new_nas[center_idx:center_idx] = gaps
    return new_nas


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def gather_gaps(
    refnas: list[NAPosition], seqnas: list[NAPosition], min_gap_distance: int
) -> tuple[list[NAPosition], list[NAPosition]]:
    """Gather nearby gaps into single windows.

    For each window, redundant gaps are removed and the remaining
    gaps are moved to the centre of the non-gap content.
    """
    slicekey: slice
    win_refnas: list[NAPosition]
    win_seqnas: list[NAPosition]
    # reverse windows so the assignment won't change index
    for slicekey in reversed(find_windows_with_gap(refnas, seqnas, min_gap_distance)):
        win_refnas = refnas[slicekey]
        win_seqnas = seqnas[slicekey]
        win_refnas, win_seqnas = remove_redundant_gaps(win_refnas, win_seqnas)

        win_refnas = move_gaps_to_center(win_refnas)
        win_seqnas = move_gaps_to_center(win_seqnas)

        refnas[slicekey] = win_refnas
        seqnas[slicekey] = win_seqnas

    return refnas, seqnas


# -----------------------------------------------------------------------
# adjust_gap_placement
# -----------------------------------------------------------------------


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def adjust_gap_placement(
    refcodons: list[list[NAPosition]],
    seqcodons: list[list[NAPosition]],
    window_size: int,
    gap_placement_score: dict[int, dict[tuple[int, int], int]],
    is_seq_start: bool,
    is_seq_end: bool,
) -> tuple[list[list[NAPosition]], list[list[NAPosition]]]:
    """Adjust gap placement for each contiguous gap group.

    Groups consecutive gap codons, extends windows left/right by up to
    *window_size* gap-free codons, then calls
    ``paired_find_best_matches`` to optimise gap position within each
    extended window.
    """
    start: int
    end: int
    offset: int
    gap_type: int
    codonpairs: Iterable[CodonPair]
    refcd: list[NAPosition]
    seqcd: list[NAPosition]
    ext_refcds: list[list[NAPosition]]
    ext_seqcds: list[list[NAPosition]]

    trim_slice: slice = find_codon_trim_slice(seqcodons)

    gap_groups: Iterable[
        tuple[
            int,  # group key: NOGAP, REFGAP or SEQGAP
            Iterable[CodonPair],
        ]
    ] = groupby(
        list(enumerate(zip(refcodons, seqcodons, strict=False)))[trim_slice],
        codon_pairs_group_key,
    )

    for gap_type, codonpairs in gap_groups:
        if gap_type == NOGAP:
            continue

        codonpairs = list(codonpairs)
        start, end = codonpairs[0][0], codonpairs[-1][0] + 1
        refcds: list[list[NAPosition]] = []
        seqcds: list[list[NAPosition]] = []
        for _, (refcd, seqcd) in codonpairs:
            refcds.append(refcd)
            seqcds.append(seqcd)

        # extend refcds/seqcds
        ext_refcds, ext_seqcds, offset = extend_codons_until_gap(
            refcodons[max(0, start - window_size) : start],
            seqcodons[max(0, start - window_size) : start],
            LEFT,
        )
        if offset:
            refcds = ext_refcds + refcds
            seqcds = ext_seqcds + seqcds
            start -= offset

        ext_refcds, ext_seqcds, offset = extend_codons_until_gap(
            refcodons[end : end + window_size],
            seqcodons[end : end + window_size],
            RIGHT,
        )
        if offset:
            refcds = refcds + ext_refcds
            seqcds = seqcds + ext_seqcds
            end += offset

        win_refnas = list(chain(*refcds))
        win_seqnas = list(chain(*seqcds))

        win_refnas, win_seqnas = paired_find_best_matches(
            win_refnas,
            win_seqnas,
            gap_type,
            gap_placement_score,
            is_seq_start and start == trim_slice.start,
            is_seq_end and end == trim_slice.stop,
        )
        (win_refcodons, win_seqcodons) = group_by_codons(win_refnas, win_seqnas)
        refcodons[start:end] = win_refcodons
        seqcodons[start:end] = win_seqcodons

    return refcodons, seqcodons


# -----------------------------------------------------------------------
# realign_gaps
# -----------------------------------------------------------------------


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def realign_gaps(
    refnas: list[NAPosition],
    seqnas: list[NAPosition],
    min_gap_distance: int,
    window_size: int,
    gap_placement_score: dict[int, dict[tuple[int, int], int]],
    is_seq_start: bool,
    is_seq_end: bool,
) -> tuple[list[NAPosition], list[NAPosition]]:
    """Gather gaps and optimise their placement.

    This is the main internal entry point that chains gap gathering,
    codon grouping, gap placement adjustment, and codon-end gap movement.

    Args:
        refnas: Reference NA positions for the active window.
        seqnas: Sequence NA positions for the active window.
        min_gap_distance: Maximum NA distance for gap merging.
        window_size: Number of flanking codons to consider.
        gap_placement_score: Bonus/penalty map.
        is_seq_start: True if the window starts at sequence start.
        is_seq_end: True if the window ends at sequence end.

    Returns:
        Updated ``(refnas, seqnas)`` with gaps optimally placed.
    """
    refcodons: list[list[NAPosition]]
    seqcodons: list[list[NAPosition]]

    refnas, seqnas = gather_gaps(refnas, seqnas, min_gap_distance)

    refcodons, seqcodons = group_by_codons(refnas, seqnas)
    refcodons, seqcodons = adjust_gap_placement(
        refcodons, seqcodons, window_size, gap_placement_score, is_seq_start, is_seq_end
    )

    # move gaps in seqcodons to codon ends
    seqcodons = move_gap_to_codon_end(seqcodons)

    return list(chain(*refcodons)), list(chain(*seqcodons))


# -----------------------------------------------------------------------
# Main entry point — pure-Python codon_align
# -----------------------------------------------------------------------


@cython.ccall
@cython.returns(tuple)
def codon_align(
    refseq: Sequence,
    seq: Sequence,
    min_gap_distance: int,
    window_size: int,
    gap_placement_score: dict[int, dict[tuple[int, int], int]],
    ref_start: int,
    ref_end: int,
) -> RefSeqPair:
    """Perform codon-aware gap realignment on a reference/sequence pair.

    Determines the active reading-frame window from *ref_start* /
    *ref_end*, checks for gaps, and delegates to ``realign_gaps`` for
    the actual optimisation.

    Args:
        refseq: Reference ``Sequence`` object.
        seq: Target ``Sequence`` object.
        min_gap_distance: Merge gaps within this NA distance.
        window_size: Flanking codons used during gap search.
        gap_placement_score: ``{gap_type: {(pos, size): score}}`` map.
        ref_start: First reference position (1-based, inclusive).
        ref_end: Last reference position (1-based, inclusive).

    Returns:
        Updated ``(refseq, seq)`` pair with codon-aligned gaps.
    """
    refnas: list[NAPosition] = refseq.seqtext
    seqnas: list[NAPosition] = seq.seqtext

    seq_idx_start: int = 0
    seq_idx_end: int = len(seqnas)

    ref_idx_start, ref_idx_end = NAPosition.posrange2indexrange(
        refnas, ref_start, ref_end, include_boundary_gaps=True
    )

    # Determine the application boundary
    # 1) follow user-specific reference boundary (ref_start, ref_end), and
    # 2) extend codon-alignment to include ref & seq boundary gaps
    # 3) ensure ref_idx_start is at the begining of codon
    idx_start: int = ref_idx_start if ref_idx_start > seq_idx_start else seq_idx_start

    while True:
        test_idx_start = NAPosition.min_nongap_index(refnas, idx_start)
        if test_idx_start < 0:
            break
        if (refnas[test_idx_start].pos - ref_start) % 3 == 0:
            break
        idx_start = test_idx_start + 1

    idx_end: int = ref_idx_end if ref_idx_end < seq_idx_end else seq_idx_end

    if idx_start == idx_end:
        # nothing to be codon aligned
        return refseq, seq

    is_seq_start: bool = idx_start <= NAPosition.min_nongap_index(seqnas)
    is_seq_end: bool = idx_end > NAPosition.max_nongap_index(seqnas)

    # step 1: apply reading frame
    refnas = refnas[idx_start:idx_end]
    seqnas = seqnas[idx_start:idx_end]

    # D2: Fast path — no gaps at all
    if not NAPosition.any_has_gap(refnas) and not NAPosition.any_has_gap(seqnas):
        return refseq, seq

    # step 2: gather and re-align nearby gaps located in same window
    refnas, seqnas = realign_gaps(
        refnas,
        seqnas,
        min_gap_distance,
        window_size,
        gap_placement_score,
        is_seq_start,
        is_seq_end,
    )

    # step 3: save "codon aligned" refseq and seq
    refseq = refseq.push_seqtext(
        refseq.seqtext[:idx_start] + refnas + refseq.seqtext[idx_end:],
        f'codonalign({ref_start},{ref_end})',
        0,
    )
    seq = seq.push_seqtext(
        seq.seqtext[:idx_start] + seqnas + seq.seqtext[idx_end:],
        f'codonalign({ref_start},{ref_end})',
        0,
    )
    return refseq, seq


# -----------------------------------------------------------------------
# CLI helpers
# -----------------------------------------------------------------------


@cython.ccall
@cython.returns(dict)
def parse_gap_placement_score(value: str) -> dict[int, dict[tuple[int, int], int]]:
    """Parse a comma-separated gap-placement-score string.

    Each token has the form ``<pos>[/<size>](ins|del):<score>``.

    Args:
        value: Raw score string, e.g. ``"204ins:-5,2041/12del:10"``.

    Returns:
        ``{REFGAP: {(pos, size): score}, SEQGAP: {…}}``.

    Raises:
        ValueError: If any token does not match the expected pattern.
    """
    pos_size: str
    gap_type: str
    gap_score: str
    score_str: str
    scores: dict[int, dict[tuple[int, int], int]] = {REFGAP: {}, SEQGAP: {}}
    for score_str in value.split(','):
        if not score_str:
            continue
        match: re.Match | None = GAP_PLACEMENT_SCORE_PATTERN.match(score_str)
        if not match:
            raise ValueError(
                'parse_gap_placement_score() is provided with an '
                f'invalid argument value: {score_str!r}'
            )
        pos_start, pos_size, gap_type, gap_score = match.groups()
        scores[REFGAP if gap_type == 'ins' else SEQGAP][
            (int(pos_start), int(pos_size) if pos_size else 0)
        ] = int(gap_score)
    return scores


@cython.ccall
@cython.returns(dict)
def gap_placement_score_callback(
    ctx: click.Context, param: click.Option, value: tuple[str]
) -> dict[int, dict[tuple[int, int], int]]:
    """Click callback that parses ``--gap-placement-score`` values."""
    if not param.name:
        raise click.BadParameter(
            'Internal error (gap_placement_score_callback:1)', ctx, param
        )
    try:
        result: dict[int, dict[tuple[int, int], int]] = parse_gap_placement_score(
            ','.join(value)
        )
        return result
    except ValueError as exp:
        raise click.BadOptionUsage(param.name, str(exp), ctx) from exp


# -----------------------------------------------------------------------
# CLI command
# -----------------------------------------------------------------------


@cli.command('codon-alignment')
@click.option(
    '--min-gap-distance',
    type=int,
    default=30,
    help=(
        'Minimal NA gap distance of the output, gaps within the '
        'minimal distance will be gathered into a single gap'
    ),
)
@click.option(
    '--window-size',
    type=int,
    default=10,
    help=(
        'AA local window size for finding the local optimal '
        'placement (BLOSUM62) for an insertion or deletion gap: '
        'the larger the window the better the result and the slower '
        'the process'
    ),
)
@click.option(
    '--gap-placement-score',
    type=str,
    multiple=True,
    default=[],
    callback=gap_placement_score_callback,
    help=(
        'Bonus (positive number) or penalty (negative number) for gaps '
        'appear at certain NA position (relative to the WHOLE ref seq) in the '
        'ref seq (ins) or target seq (del). For example, 204ins:-5 is a -5 '
        'penalty designate to a gap with any size gap in ref seq after NA '
        'position 204 (AA position 68). 2041/12del:10 is a +10 score for a '
        '12 NAs size (4 codons) gap in target seq at NA position 2041, '
        'equivalent to deletion at 681, 682, 683 and 684 AA position. '
        'Multiple scores can be delimited by commas, such as '
        '204ins:-5,2041/12del:10.'
    ),
)
@click.option(
    '--backend',
    type=click.Choice(['python', 'rust'], case_sensitive=False),
    default='rust',
    help=(
        'Execution backend: "python" uses the pure-Python/Cython '
        'implementation; "rust" (default) uses the Rust-accelerated '
        'implementation via postalign_rs.'
    ),
)
@click.argument('ref_start', type=int, default=1)
@click.argument('ref_end', type=int, default=-1)
def codon_alignment(
    min_gap_distance: int,
    window_size: int,
    # For NASize, 0 means any size
    #                        Indel         NAPos NASize Score
    #                          v              v     v     v
    gap_placement_score: dict[int, dict[tuple[int, int], int]],
    backend: str,
    ref_start: int,
    ref_end: int,
    # XXX: see https://github.com/cython/cython/issues/2753
    # this has been fixed by cython 3.0
    # ) -> Processor[Iterable[RefSeqPair]]:
) -> Processor:
    """Perform codon alignment.

    A re-implementation of the "codon-align" tool created by LANL HIV
    Sequence Database.  Supports a pure-Python/Cython backend and a
    Rust-accelerated backend selectable via ``--backend``.

    The arguments <REF_START> and <REF_END> provide the position range
    (relative to ref. sequence) where the codon alignment should be applied.
    """
    if ref_start < 1:
        raise click.ClickException(
            f'argument <REF_START>:{ref_start} must be not less than 1'
        )
    if ref_end > 0 and ref_end - 2 < ref_start:
        raise click.ClickException(
            f'no enough codon between arguments <REF_START>:{ref_start} and <REF_END>:{ref_end}'
        )

    # Select backend implementation
    if backend == 'rust':
        from .codon_alignment_rust import codon_align as _align_fn
    else:
        _align_fn = codon_align

    @intermediate_processor('codon-alignment')
    def processor(iterator: Iterable[RefSeqPair], *args: Any) -> Iterable[RefSeqPair]:
        refseq: Sequence
        seq: Sequence
        for refseq, seq in iterator:
            if refseq.seqtype != NAPosition:
                raise click.ClickException(
                    'Codon alignment only applies to nucleotide sequences.'
                )

            seqnas: list[NAPosition] = refseq.seqtext

            if not seqnas:
                # skip empty sequences
                yield refseq, seq
            else:
                my_ref_end: int = ref_end
                if my_ref_end <= 0:
                    my_ref_end = NAPosition.max_pos(seqnas)
                yield _align_fn(
                    refseq,
                    seq,
                    min_gap_distance,
                    window_size,
                    gap_placement_score,
                    ref_start,
                    my_ref_end,
                )

    return processor
