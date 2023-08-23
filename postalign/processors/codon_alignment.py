import re
import click
import cython  # type: ignore
from typing import Iterable, Tuple, List, Set, Optional, Dict, Any
from itertools import chain, groupby

from ..cli import cli
from ..utils import group_by_codons, find_codon_trim_slice
from ..models import Sequence, RefSeqPair, NAPosition
from ..utils.codonutils import translate_codons
from ..utils.iupac import iupac_score
from ..utils.blosum62 import blosum62_score

from ..processor import intermediate_processor, Processor


NOGAP: int = 0b00
REFGAP: int = 0b01
SEQGAP: int = 0b10

LEFT: int = 0b00
RIGHT: int = 0b01

GAP_PLACEMENT_SCORE_PATTERN: re.Pattern = re.compile(
    r'^(\d+)(?:/(\d+))?(ins|del):(-?\d+)$'
)

CodonPair = Tuple[
    int,  # refpos0
    Tuple[
        List[NAPosition],  # refcodon
        List[NAPosition]   # poscodon
    ]
]


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def extend_codons_until_gap(
    ref_codons: List[List[NAPosition]],
    seq_codons: List[List[NAPosition]],
    direction: int
) -> Tuple[
    List[List[NAPosition]],
    List[List[NAPosition]],
    int
]:
    if direction == LEFT:
        ref_codons.reverse()
        seq_codons.reverse()
    idx: int
    refcd: List[NAPosition]
    seqcd: List[NAPosition]
    endidx: int = len(ref_codons)
    for idx, (refcd, seqcd) in enumerate(zip(ref_codons, seq_codons)):
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
    refnas: List[NAPosition],
    seqnas: List[NAPosition],
    min_gap_distance: int
) -> List[slice]:
    """List all windows with gap:

    Examples:

    ref: AAAA-----------
    seq: -----------BBBB

    ref: -----------
    seq: -----------

    ref: ---AAAA----
    seq: ---BBBB----
    """
    refna: NAPosition
    seqna: NAPosition
    first_gap_idx: int = -1
    last_gap_idx: int = -1
    windows: List[slice] = []
    for idx, (refna, seqna) in enumerate(zip(refnas, seqnas)):
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
def find_first_gap(nas: List[NAPosition]) -> int:
    idx: int
    na: NAPosition
    for idx, na in enumerate(nas):
        if na.is_gap:
            return idx
    return -1


@cython.cfunc
@cython.inline
@cython.returns(list)
def move_gap_to_codon_end(
    codons: List[List[NAPosition]]
) -> List[List[NAPosition]]:
    new_codons: List[List[NAPosition]] = []
    for codon in codons:
        new_codons.append(
            [na for na in codon if not na.is_gap] +
            [na for na in codon if na.is_gap]
        )
    return new_codons


@cython.cfunc
@cython.inline
def calc_match_score(
    mynas: List[NAPosition],
    othernas: List[NAPosition],
    base_score: float
) -> float:
    myna: NAPosition
    otherna: NAPosition
    myaa: bytes
    otheraa: bytes
    myaas: List[bytes] = translate_codons(mynas)
    otheraas: List[bytes] = translate_codons(othernas)
    score: float = base_score
    for myna, otherna in zip(mynas, othernas):
        score += iupac_score(myna.notation, otherna.notation)
    for myaa, otheraa in zip(myaas, otheraas):
        score += blosum62_score(myaa, otheraa)
    return score


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def separate_gaps_from_nas(
    nas: List[NAPosition]
) -> Tuple[List[NAPosition], List[NAPosition]]:
    na: NAPosition
    nongaps: List[NAPosition] = []
    gaps: List[NAPosition] = []
    for na in nas:
        if na.is_gap:
            gaps.append(na)
        else:
            nongaps.append(na)
    return nongaps, gaps


@cython.cfunc
@cython.inline
@cython.returns(list)
def find_best_matches(
    mynas: List[NAPosition],
    othernas: List[NAPosition],
    bp1_indices: Set[int],
    gap_type: int,
    gap_placement_score: Dict[Tuple[int, int], int],
    is_start: bool,
    is_end: bool
) -> List[NAPosition]:
    idx: int
    score: Tuple[float, int]
    mygap: List[NAPosition]
    test_mynas: List[NAPosition]
    orig_gapidx: int = find_first_gap(mynas)
    mynas, mygap = separate_gaps_from_nas(mynas)
    gaplen: int = len(mygap)
    max_score: Optional[Tuple[float, int]] = None
    best_mynas: Optional[List[NAPosition]] = None
    scanstart: int = 3 if gap_type == REFGAP else 0
    mynas_len: int = len(mynas)
    for idx in range(scanstart, mynas_len + 1, 3):
        napos: int
        test_mynas = mynas[::]
        test_mynas[idx:idx] = mygap
        base_score: float = float(-gaplen)
        if is_start and idx == 0:
            base_score = .0
        elif is_end and idx + 3 > mynas_len:
            base_score = .0
        score_val: float = calc_match_score(test_mynas, othernas, base_score)
        if gap_type == REFGAP:
            napos = mynas[idx - 1].pos
        else:  # gap_type == SEQGAP
            napos = othernas[idx].pos
        if (napos, gaplen) in gap_placement_score:
            score_val += gap_placement_score[(napos, gaplen)]
        elif (napos, 0) in gap_placement_score:
            score_val += gap_placement_score[(napos, 0)]

        if idx in bp1_indices:
            # reward gaps inserted between codons
            score = (score_val + 1, 2)
        elif idx == orig_gapidx:
            # respect the original gapidx if it's already one of the best
            score = (score_val, 1)
        else:
            score = (score_val, 0)
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
    refnas: List[NAPosition],
    seqnas: List[NAPosition],
    gap_type: int,
    gap_placement_score: Dict[int, Dict[Tuple[int, int], int]],
    is_seq_start: bool,
    is_seq_end: bool
) -> Tuple[List[NAPosition], List[NAPosition]]:
    idx: int
    na: NAPosition
    bp1_indices: Set[int] = set()
    bp: int = 0
    for idx, na in enumerate(refnas):
        if na.is_gap:
            continue
        bp = (bp + 1) % 3
        if bp == 1:
            bp1_indices.add(idx)

    if gap_type == REFGAP:
        refnas = find_best_matches(
            refnas, seqnas, bp1_indices,
            gap_type,
            gap_placement_score[gap_type],
            # for REFGAPs, ending gaps also have penalty
            False,
            False)
    elif gap_type == SEQGAP:
        seqnas = find_best_matches(
            seqnas, refnas, bp1_indices,
            gap_type,
            gap_placement_score[gap_type],
            is_seq_start,
            is_seq_end)
    return refnas, seqnas


@cython.ccall
def codon_pairs_group_key(cdpair: Tuple[
    int,
    Tuple[List[NAPosition], List[NAPosition]]
]) -> int:
    refcd: List[NAPosition]
    seqcd: List[NAPosition]
    _, (refcd, seqcd) = cdpair
    if NAPosition.any_has_gap(refcd):
        return REFGAP
    if NAPosition.any_has_gap(seqcd):
        return SEQGAP
    return NOGAP


@cython.cfunc
@cython.inline
@cython.returns(list)
def move_gaps_to_center(nas: List[NAPosition]) -> List[NAPosition]:
    new_nas: List[NAPosition]
    gaps: List[NAPosition]
    new_nas, gaps = separate_gaps_from_nas(nas)
    center_idx: int = len(new_nas) // 2
    new_nas[center_idx:center_idx] = gaps
    return new_nas


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def gather_gaps(
    refnas: List[NAPosition],
    seqnas: List[NAPosition],
    min_gap_distance: int
) -> Tuple[
    List[NAPosition],
    List[NAPosition]
]:
    """Gather gaps together according to window"""
    slicekey: slice
    win_refnas: List[NAPosition]
    win_seqnas: List[NAPosition]
    # reverse windows so the assignment won't change index
    for slicekey in reversed(find_windows_with_gap(
        refnas, seqnas, min_gap_distance
    )):
        win_refnas = refnas[slicekey]
        win_seqnas = seqnas[slicekey]
        win_refnas, win_seqnas = remove_redundant_gaps(win_refnas, win_seqnas)

        win_refnas = move_gaps_to_center(win_refnas)
        win_seqnas = move_gaps_to_center(win_seqnas)

        refnas[slicekey] = win_refnas
        seqnas[slicekey] = win_seqnas

    return refnas, seqnas


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def adjust_gap_placement(
    refcodons: List[List[NAPosition]],
    seqcodons: List[List[NAPosition]],
    window_size: int,
    gap_placement_score: Dict[int, Dict[Tuple[int, int], int]],
    is_seq_start: bool,
    is_seq_end: bool
) -> Tuple[List[List[NAPosition]], List[List[NAPosition]]]:
    """Adjust continuous refgap/seqgap placement"""
    start: int
    end: int
    offset: int
    gap_type: int
    codonpairs: Iterable[CodonPair]
    refcd: List[NAPosition]
    seqcd: List[NAPosition]
    ext_refcds: List[List[NAPosition]]
    ext_seqcds: List[List[NAPosition]]

    trim_slice: slice = find_codon_trim_slice(seqcodons)

    gap_groups: Iterable[
        Tuple[
            int,  # group key: NOGAP, REFGAP or SEQGAP
            Iterable[CodonPair]
        ]
    ] = groupby(
        list(enumerate(zip(refcodons, seqcodons)))[trim_slice],
        codon_pairs_group_key
    )

    for gap_type, codonpairs in gap_groups:
        if gap_type == NOGAP:
            continue

        codonpairs = list(codonpairs)
        start, end = codonpairs[0][0], codonpairs[-1][0] + 1
        refcds: List[List[NAPosition]] = []
        seqcds: List[List[NAPosition]] = []
        for _, (refcd, seqcd) in codonpairs:
            refcds.append(refcd)
            seqcds.append(seqcd)

        # extend refcds/seqcds
        ext_refcds, ext_seqcds, offset = extend_codons_until_gap(
            refcodons[max(0, start - window_size):start],
            seqcodons[max(0, start - window_size):start],
            LEFT
        )
        if offset:
            refcds = ext_refcds + refcds
            seqcds = ext_seqcds + seqcds
            start -= offset

        ext_refcds, ext_seqcds, offset = extend_codons_until_gap(
            refcodons[end:end + window_size],
            seqcodons[end:end + window_size],
            RIGHT
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
            is_seq_end and end == trim_slice.stop
        )
        (win_refcodons,
         win_seqcodons) = group_by_codons(win_refnas, win_seqnas)
        refcodons[start:end] = win_refcodons
        seqcodons[start:end] = win_seqcodons

    return refcodons, seqcodons


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def realign_gaps(
    refnas: List[NAPosition],
    seqnas: List[NAPosition],
    min_gap_distance: int,
    window_size: int,
    gap_placement_score: Dict[int, Dict[Tuple[int, int], int]],
    is_seq_start: bool,
    is_seq_end: bool
) -> Tuple[List[NAPosition], List[NAPosition]]:
    refcodons: List[List[NAPosition]]
    seqcodons: List[List[NAPosition]]

    refnas, seqnas = gather_gaps(refnas, seqnas, min_gap_distance)

    refcodons, seqcodons = group_by_codons(refnas, seqnas)
    refcodons, seqcodons = adjust_gap_placement(
        refcodons,
        seqcodons,
        window_size,
        gap_placement_score,
        is_seq_start,
        is_seq_end
    )

    # move gaps in seqcodons to codon ends
    seqcodons = move_gap_to_codon_end(seqcodons)

    return list(chain(*refcodons)), list(chain(*seqcodons))


@cython.cfunc
@cython.inline
@cython.returns(list)
def remove_n_gaps(nas: List[NAPosition], n_gaps: int) -> List[NAPosition]:
    na: NAPosition
    result: List[NAPosition] = []
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
    refnas: List[NAPosition],
    seqnas: List[NAPosition]
) -> Tuple[List[NAPosition], List[NAPosition]]:
    """Remove redundant gaps

    Gap should only exist in either sequence but not both

    Example of matched gaps:

                vv
    TCAGCTGATGCA---CAA
    TCAGGC-TGATGCAC-AA
          ^        ^

    Two of above gaps are redundant and should be removed
    before the alignment get further processed
    """
    n_gaps: int = min(
        NAPosition.count_gaps(refnas),
        NAPosition.count_gaps(seqnas)
    )
    if n_gaps:
        refnas = remove_n_gaps(refnas, n_gaps)
        seqnas = remove_n_gaps(seqnas, n_gaps)
    return refnas, seqnas


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
    # 1) follow user-specific reference boundary (ref_start, ref_end), and
    # 2) extend codon-alignment to include ref & seq boundary gaps
    # 3) ensure ref_idx_start is at the begining of codon
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

    # step 2: gather and re-align nearby gaps located in same window
    refnas, seqnas = realign_gaps(
        refnas,
        seqnas,
        min_gap_distance,
        window_size,
        gap_placement_score,
        is_seq_start,
        is_seq_end)

    # step 3: save "codon aligned" refseq and seq
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


@cython.ccall
@cython.returns(dict)
def parse_gap_placement_score(value: str) -> Dict[
    int, Dict[Tuple[int, int], int]
]:
    pos: str
    pos_size: str
    gap_type: str
    gap_score: str
    score_str: str
    scores: Dict[int, Dict[Tuple[int, int], int]] = {
        REFGAP: {},
        SEQGAP: {}
    }
    for score_str in value.split(','):
        if not score_str:
            continue
        match: Optional[re.Match] = \
            GAP_PLACEMENT_SCORE_PATTERN.match(score_str)
        if not match:
            raise ValueError(
                'parse_gap_placement_score() is provided with an '
                'invalid argument value: {!r}'.format(score_str)
            )
        pos_start, pos_size, gap_type, gap_score = match.groups()
        scores[REFGAP if gap_type == 'ins' else SEQGAP][(
            int(pos_start),
            int(pos_size) if pos_size else 0
        )] = int(gap_score)
    return scores


@cython.ccall
@cython.returns(dict)
def gap_placement_score_callback(
    ctx: click.Context,
    param: click.Option,
    value: Tuple[str]
) -> Dict[int, Dict[Tuple[int, int], int]]:
    if not param.name:
        raise click.BadParameter(
            'Internal error (gap_placement_score_callback:1)',
            ctx,
            param
        )
    try:
        result: Dict[
            int, Dict[Tuple[int, int], int]
        ] = parse_gap_placement_score(','.join(value))
        return result
    except ValueError as exp:
        raise click.BadOptionUsage(
            param.name,
            str(exp),
            ctx
        )


@cli.command('codon-alignment')
@click.option(
    '--min-gap-distance',
    type=int,
    default=30,
    help=(
        'Minimal NA gap distance of the output, gaps within the'
        'minnimal distance will be gathered into a single gap'
    ))
@click.option(
    '--window-size',
    type=int,
    default=10,
    help=(
        'AA local window size for finding the local optimal '
        'placement (BLOSUM62) for an insertion or deletion gap: '
        'the larger the window the better the result and the slower '
        'the process'
    ))
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
    ))
@click.argument(
    'ref_start', type=int, default=1
)
@click.argument(
    'ref_end', type=int, default=-1
)
def codon_alignment(
    min_gap_distance: int,
    window_size: int,
    # For NASize, 0 means any size
    #                        Indel         NAPos NASize Score
    #                          v              v     v     v
    gap_placement_score: Dict[int, Dict[Tuple[int, int], int]],
    ref_start: int,
    ref_end: int
    # XXX: see https://github.com/cython/cython/issues/2753
    # this has been fixed by cython 3.0
    # ) -> Processor[Iterable[RefSeqPair]]:
) -> Processor:
    """Codon alignment

    A blackbox re-implementation of the "codon-align" tool
    created by LANL HIV Sequence Database.

    The arguments <REF_START> and <REF_END> provide the position range
    (relative to ref. sequence) where the codon alignment should be applied.
    """
    if ref_start < 1:
        raise click.ClickException(
            'argument <REF_START>:{} must be not less than 1'.format(ref_start)
        )
    if ref_end > 0 and ref_end - 2 < ref_start:
        raise click.ClickException(
            'no enough codon between arguments <REF_START>:{} and <REF_END>:{}'
            .format(ref_start, ref_end)
        )

    @intermediate_processor('codon-alignment')
    def processor(
        iterator: Iterable[RefSeqPair], *args: Any
    ) -> Iterable[RefSeqPair]:
        refseq: Sequence
        seq: Sequence
        for refseq, seq in iterator:
            if refseq.seqtype != NAPosition:
                raise click.ClickException(
                    'Codon alignment only applies to nucleotide '
                    'sequences.')

            seqnas: List[NAPosition] = refseq.seqtext

            if not seqnas:
                # skip empty sequences
                yield refseq, seq
            else:
                my_ref_end: int = ref_end
                if my_ref_end <= 0:
                    my_ref_end = NAPosition.max_pos(seqnas)
                yield codon_align(
                    refseq,
                    seq,
                    min_gap_distance,
                    window_size,
                    gap_placement_score,
                    ref_start, my_ref_end
                )

    return processor
