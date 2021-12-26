import click
import cython  # type: ignore
from typing import Iterable, Tuple, List, Set, Optional
from itertools import chain, groupby

from ..cli import cli
from ..utils import group_by_codons
from ..models import Sequence, RefSeqPair, NAPosition, NAPosOrList
from ..utils.codonutils import translate_codons
from ..utils.blosum62 import blosum62_score

from ..processor import intermediate_processor, Processor


NOGAP: int = 0b00
REFGAP: int = 0b01
SEQGAP: int = 0b10

LEFT: int = 0b00
RIGHT: int = 0b01

CodonPair = Tuple[
    int,  # refpos0
    Tuple[
        List[NAPosition],  # refcodon
        List[NAPosition]   # poscodon
    ]
]


@cython.cfunc
@cython.inline
def find_first_gap(
    nas: List[NAPosition],
    start: int = 0
) -> int:
    idx: int
    for idx, na in enumerate(nas[start:]):
        if na.is_single_gap:
            return start + idx
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
            [na for na in codon if not na.is_single_gap] +
            [na for na in codon if na.is_single_gap]
        )
    return new_codons


@cython.cfunc
@cython.inline
def calc_match_score(
    mynas: List[NAPosition],
    othernas: List[NAPosition]
) -> float:
    myaa: bytes
    otheraa: bytes
    myaas: List[bytes] = translate_codons(mynas)
    otheraas: List[bytes] = translate_codons(othernas)
    score: float = .0
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
    nongaps: List[NAPosition] = [na for na in nas if not na.is_single_gap]
    gaps: List[NAPosition] = [na for na in nas if na.is_single_gap]
    return nongaps, gaps


@cython.cfunc
@cython.inline
@cython.returns(list)
def find_best_matches(
    mynas: List[NAPosition],
    othernas: List[NAPosition],
    bp1_indices: Set[int],
    scanstart: int = 0,
    scanstep: int = 1
) -> List[NAPosition]:
    idx: int
    score: Tuple[float, int]
    mygap: List[NAPosition]
    test_mynas: List[NAPosition]
    orig_gapidx: int = find_first_gap(mynas)
    mynas, mygap = separate_gaps_from_nas(mynas)
    max_score: Optional[Tuple[float, int]] = None
    best_mynas: Optional[List[NAPosition]] = None
    for idx in range(scanstart, len(mynas) + 1, scanstep):
        test_mynas = mynas[::]
        test_mynas[idx:idx] = mygap
        score_val: float = calc_match_score(test_mynas, othernas)
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
    gap_type: int
) -> Tuple[List[NAPosition], List[NAPosition]]:
    idx: int
    na: NAPosition
    bp1_indices: Set[int] = set()
    bp: int = 0
    for idx, na in enumerate(refnas):
        if na.is_single_gap:
            continue
        bp = (bp + 1) % 3
        if bp == 1:
            bp1_indices.add(idx)

    if gap_type == REFGAP:
        refnas = find_best_matches(refnas, seqnas, bp1_indices,
                                   scanstart=3, scanstep=3)
    elif gap_type == SEQGAP:
        seqnas = find_best_matches(seqnas, refnas, bp1_indices,
                                   scanstart=0, scanstep=3)
    return refnas, seqnas


@cython.ccall
def codon_pairs_group_key(cdpair: Tuple[
    int,
    Tuple[List[NAPosition], List[NAPosition]]
]) -> int:
    refcd: List[NAPosition]
    seqcd: List[NAPosition]
    _, (refcd, seqcd) = cdpair
    if NAPosition.list_contains_any_gap(refcd):
        return REFGAP
    if NAPosition.list_contains_any_gap(seqcd):
        return SEQGAP
    return NOGAP


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
            if na.is_single_gap:
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
def move_gaps_to_index(nas: List[NAPosition], index: int) -> List[NAPosition]:
    new_nas: List[NAPosition]
    gaps: List[NAPosition]
    new_nas, gaps = separate_gaps_from_nas(nas)
    new_nas[index:index] = gaps
    return new_nas


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def gather_gaps(
    refcodons: List[List[NAPosition]],
    seqcodons: List[List[NAPosition]],
    window_size: int
) -> Tuple[
    List[List[NAPosition]],
    List[List[NAPosition]]
]:
    """Gather gaps together according to window"""
    pointer: int
    slicekey: slice
    win_refcodons: List[List[NAPosition]]
    win_seqcodons: List[List[NAPosition]]
    win_refnas: List[NAPosition]
    win_seqnas: List[NAPosition]
    pointer_limit: int = len(refcodons) - window_size
    for pointer in range(0, pointer_limit):
        slicekey = slice(pointer, pointer + window_size)
        win_refcodons = refcodons[slicekey]
        win_seqcodons = seqcodons[slicekey]
        win_refnas = list(chain(*win_refcodons))
        win_seqnas = list(chain(*win_seqcodons))
        gapidx: int = find_first_gap(win_refnas)

        if gapidx > -1:
            # When gap(s) are found in refnas, first move all gaps in a window
            # to the first gap position, for example:
            # ref: AAA-CCCT--TT  => AAA---CCCTTT
            win_refnas = move_gaps_to_index(win_refnas, gapidx)

            # Then align gaps in seqnas with refnas, for example:
            # ref: AAA---CCCTTT     AAA---CCCTTT
            #                    =>
            # seq: AAACCC---TTT     AAA---CCCTTT
            win_seqnas = move_gaps_to_index(win_seqnas, gapidx)
        else:
            gapidx = find_first_gap(win_seqnas)
            if gapidx > -1:
                # When gap(s) are found in seqnas, move all gaps in a window
                # to the first gap position like above did for refnas
                win_seqnas = move_gaps_to_index(win_seqnas, gapidx)

        win_refnas, win_seqnas = remove_matching_gaps(win_refnas, win_seqnas)
        (win_refcodons,
         win_seqcodons) = group_by_codons(win_refnas, win_seqnas)
        refcodons[slicekey] = win_refcodons
        seqcodons[slicekey] = win_seqcodons

    return refcodons, seqcodons


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def adjust_gap_placement(
    refcodons: List[List[NAPosition]],
    seqcodons: List[List[NAPosition]],
    window_size: int
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

    gap_groups: Iterable[
        Tuple[
            int,  # group key: NOGAP, REFGAP or SEQGAP
            Iterable[CodonPair]
        ]
    ] = groupby(
        enumerate(zip(refcodons,
                      seqcodons)),
        codon_pairs_group_key)

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
        win_refnas, win_seqnas = \
            paired_find_best_matches(win_refnas, win_seqnas, gap_type)
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
    window_size: int
) -> Tuple[List[NAPosition], List[NAPosition]]:
    refcodons: List[List[NAPosition]]
    seqcodons: List[List[NAPosition]]

    refcodons, seqcodons = group_by_codons(refnas, seqnas)

    refcodons, seqcodons = gather_gaps(refcodons, seqcodons, window_size)
    refcodons, seqcodons = adjust_gap_placement(
        refcodons, seqcodons, window_size)

    # move gaps in seqcodons to codon ends
    seqcodons = move_gap_to_codon_end(seqcodons)

    return list(chain(*refcodons)), list(chain(*seqcodons))


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def remove_matching_gaps(
    refnas: NAPosOrList,
    seqnas: NAPosOrList
) -> Tuple[List[NAPosition], List[NAPosition]]:
    """Remove matching gaps

    Example of matched gaps:

    TCAGCTGATGCA---CAA
    TCAGCTGATGCAC---AA
                 ^^
                 These two are "matched gaps" and can be
                 removed from pairwise alignment
    """
    refna: NAPosition
    seqna: NAPosition
    cleaned_refnas: List[NAPosition] = []
    cleaned_seqnas: List[NAPosition] = []
    for refna, seqna in zip(refnas, seqnas):
        if not refna.is_single_gap or not seqna.is_single_gap:
            cleaned_refnas.append(refna)
            cleaned_seqnas.append(seqna)
    return cleaned_refnas, cleaned_seqnas


@cython.ccall
@cython.returns(tuple)
def codon_align(
    refseq: Sequence,
    seq: Sequence,
    window_size: int,
    refstart: int,
    refend: int,
    check_boundary: bool = True
) -> RefSeqPair:
    reftext: NAPosition = refseq.seqtext
    if check_boundary and (
        reftext[0].is_single_gap or
        reftext[-1].is_single_gap
    ):
        raise click.ClickException(
            'Unable to perform codon-alignment without the alignments '
            'being trimmed properly. Can be solved by pre-processing '
            'the alignments by command "trim-by-ref".'
        )

    idxstart: int
    idxend: int
    idxstart, idxend = reftext.posrange2indexrange(refstart, refend)
    if idxstart == idxend:
        # nothing to be codon aligned
        return refseq, seq

    # step 1: apply reading frame
    reftext = reftext[idxstart:idxend]
    seqtext = seq.seqtext[idxstart:idxend]

    # step 2: remove matched gaps (due to MSA) from ref and seq
    refnas: List[NAPosition]
    seqnas: List[NAPosition]
    refnas, seqnas = remove_matching_gaps(reftext, seqtext)

    # step 3: gather and re-align nearby gaps located in same window
    refnas, seqnas = realign_gaps(refnas, seqnas, window_size)

    # last step: save "codon aligned" refseq and seq
    refseq = refseq.push_seqtext(
        refseq.seqtext[:idxstart] +
        NAPosition.join(refnas) +
        refseq.seqtext[idxend:],
        'codonalign({},{})'.format(refstart, refend), 0)
    seq = seq.push_seqtext(
        seq.seqtext[:idxstart] +
        NAPosition.join(seqnas) +
        seq.seqtext[idxend:],
        'codonalign({},{})'.format(refstart, refend), 0)
    return refseq, seq


@cli.command('codon-alignment')
@click.option(
    '--window-size',
    type=int,
    default=5,
    help=(
        'Window size as # of codons for frameshift compensation'
    ))
@click.argument(
    'ref_start', type=int, default=1
)
@click.argument(
    'ref_end', type=int, default=-1
)
def codon_alignment(
    window_size: int,
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
    def processor(iterator: Iterable[RefSeqPair]) -> Iterable[RefSeqPair]:
        refseq: Sequence
        seq: Sequence
        for refseq, seq in iterator:
            if refseq.seqtype != 'NA':
                raise click.ClickException(
                    'Codon alignment only applies to nucleotide '
                    'sequences.')

            seqtext: NAPosition = refseq.seqtext

            if seqtext.empty():
                # skip empty sequences
                yield refseq, seq
            else:
                my_ref_end: int = ref_end
                if my_ref_end <= 0:
                    my_ref_end = seqtext.max_pos()
                yield codon_align(
                    refseq,
                    seq,
                    window_size,
                    ref_start, my_ref_end
                )

    return processor
