"""Rust-accelerated codon alignment (T4).

Same codon_align() API as codon_alignment.py but the core realign_gaps
pipeline is executed in Rust via PyO3 (postalign_rs).  Includes Option D
improvements (D1 precomputed AAs, D4 center-expand search) baked into the
Rust code.
"""
import cython  # type: ignore
from typing import Tuple, List, Dict

from ..models import Sequence, RefSeqPair, NAPosition

import postalign_rs  # Rust extension


# -----------------------------------------------------------------------
# parse_gap_placement_score — reused from the Python implementation
# -----------------------------------------------------------------------

REFGAP: int = 0b01
SEQGAP: int = 0b10


def parse_gap_placement_score(
    config: str
) -> Dict[int, Dict[Tuple[int, int], int]]:
    """Parse gap_placement_score config string.

    Format: <gap_type>:<pos>:<gaplen>:<score>[,...]
    gap_type: 'R' for REFGAP, 'S' for SEQGAP
    """
    result: Dict[int, Dict[Tuple[int, int], int]] = {
        REFGAP: {},
        SEQGAP: {},
    }
    if not config:
        return result
    for item in config.split(','):
        parts = item.strip().split(':')
        if len(parts) != 4:
            continue
        gap_type_str, pos_str, gaplen_str, score_str = parts
        gap_type = REFGAP if gap_type_str == 'R' else SEQGAP
        result[gap_type][(int(pos_str), int(gaplen_str))] = int(score_str)
    return result


# -----------------------------------------------------------------------
# Main entry point — same API as original
# -----------------------------------------------------------------------

@cython.ccall
@cython.returns(tuple)
def codon_align(
    refseq: Sequence,
    seq: Sequence,
    min_gap_distance: int,
    window_size: int,
    gap_placement_score: Dict[int, Dict[Tuple[int, int], int]],
    ref_start: int,
    ref_end: int
) -> RefSeqPair:
    refnas: List[NAPosition] = refseq.seqtext
    seqnas: List[NAPosition] = seq.seqtext

    seq_idx_start: int = 0
    seq_idx_end: int = len(seqnas)

    ref_idx_start, ref_idx_end = NAPosition.posrange2indexrange(
        refnas, ref_start, ref_end, include_boundary_gaps=True)

    # Determine the application boundary
    idx_start: int = (
        ref_idx_start
        if ref_idx_start > seq_idx_start
        else seq_idx_start
    )

    while True:
        test_idx_start = NAPosition.min_nongap_index(refnas, idx_start)
        if test_idx_start < 0:
            break
        if (refnas[test_idx_start].pos - ref_start) % 3 == 0:
            break
        idx_start = test_idx_start + 1

    idx_end: int = (
        ref_idx_end
        if ref_idx_end < seq_idx_end
        else seq_idx_end
    )

    if idx_start == idx_end:
        # nothing to be codon aligned
        return refseq, seq

    is_seq_start: bool = idx_start <= NAPosition.min_nongap_index(seqnas)
    is_seq_end: bool = idx_end > NAPosition.max_nongap_index(seqnas)

    # step 1: apply reading frame
    refnas = refnas[idx_start:idx_end]
    seqnas = seqnas[idx_start:idx_end]

    if not NAPosition.any_has_gap(refnas) and \
            not NAPosition.any_has_gap(seqnas):
        return refseq, seq

    # step 2: convert NAPosition lists → flat arrays for Rust
    ref_notations = [na.notation for na in refnas]
    ref_positions = [na.pos for na in refnas]
    ref_flags = [na.flag for na in refnas]
    seq_notations = [na.notation for na in seqnas]
    seq_positions = [na.pos for na in seqnas]
    seq_flags = [na.flag for na in seqnas]

    # step 3: call Rust realign_gaps
    (out_ref_n, out_ref_p, out_ref_f,
     out_seq_n, out_seq_p, out_seq_f) = postalign_rs.realign_gaps(
        ref_notations, ref_positions, ref_flags,
        seq_notations, seq_positions, seq_flags,
        min_gap_distance,
        window_size,
        gap_placement_score,
        is_seq_start,
        is_seq_end,
    )

    # step 4: reconstruct NAPosition lists from Rust output
    refnas = [
        NAPosition(n, p, f)
        for n, p, f in zip(out_ref_n, out_ref_p, out_ref_f)
    ]
    seqnas = [
        NAPosition(n, p, f)
        for n, p, f in zip(out_seq_n, out_seq_p, out_seq_f)
    ]

    # step 5: save "codon aligned" refseq and seq
    refseq = refseq.push_seqtext(
        refseq.seqtext[:idx_start] +
        refnas +
        refseq.seqtext[idx_end:],
        'codonalign({},{})'.format(ref_start, ref_end), 0)
    seq = seq.push_seqtext(
        seq.seqtext[:idx_start] +
        seqnas +
        seq.seqtext[idx_end:],
        'codonalign({},{})'.format(ref_start, ref_end), 0)
    return refseq, seq
