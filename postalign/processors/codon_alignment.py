import click
from typing import Iterable, List, Any, Callable, Optional
from itertools import chain, groupby

from ..cli import cli
from ..utils import group_by_codons
from ..models.sequence import RefSeqPair, PositionalSeqStr
from ..utils.codonutils import translate_codons
from ..utils.blosum62 import blosum62_score

from .base import intermediate_processor, Processor


def list_index(
    listdata: List[Any],
    cond: Callable[[Any], bool],
    start: int = 0
) -> int:
    idx: int
    for idx, n in enumerate(listdata[start:]):
        # TODO: heave function calls, use inline cython
        if cond(n):
            return start + idx
    raise ValueError(
        'given condition {!r} is not found in list'
        .format(cond)
    )


def list_rindex(
    listdata: List[Any],
    cond: Callable[[Any], bool],
    start: Optional[int] = None
) -> int:
    if start is None:
        start = len(listdata)
    idx: int
    for idx, n in enumerate(listdata[:start][::-1]):
        # TODO: heave function calls, use inline cython
        if cond(n):
            return start - idx - 1
    raise ValueError(
        'given condition {!r} is not found in list'
        .format(cond)
    )


def move_gap_to_codon_end(
    codons: List[List[PositionalSeqStr]]
) -> List[List[PositionalSeqStr]]:
    new_codons: List[List[PositionalSeqStr]] = []
    for codon in codons:
        new_codons.append(
            [na for na in codon if not na.is_gap()] +
            [na for na in codon if na.is_gap()]
        )
    return new_codons


def calc_match_score(
    mynas: List[PositionalSeqStr],
    othernas: List[PositionalSeqStr]
) -> int:
    myaa: str
    otheraa: str
    myaas: str = translate_codons(mynas)
    otheraas: str = translate_codons(othernas)
    score: int = 0
    for myaa, otheraa in zip(myaas, otheraas):
        score += blosum62_score(myaa, otheraa)
    return score


def find_best_matches(mynas, othernas, bp1_indices, scanstart=0, scanstep=1):
    orig_gapidx = list_index(mynas, lambda na: na.is_gap())
    mygap = [na for na in mynas if na.is_gap()]
    mynas = [na for na in mynas if not na.is_gap()]
    max_score = None
    best_mynas = None
    for idx in range(scanstart, len(mynas) + 1, scanstep):
        test_mynas = mynas[::]
        test_mynas[idx:idx] = mygap
        score = calc_match_score(test_mynas, othernas)
        if idx in bp1_indices:
            # reward gaps inserted between codons
            score = (score + 1, 2)
        elif idx == orig_gapidx:
            # respect the original gapidx if it's already one of the best
            score = (score, 1)
        else:
            score = (score, 0)
        if max_score is None or score > max_score:
            max_score = score
            best_mynas = test_mynas
    if best_mynas is None:
        # fallback to mynas, if no best match is found
        best_mynas = mynas
    return best_mynas


def paired_find_best_matches(refnas, seqnas):
    bp1_indices = set()
    bp = 0
    for idx, na in enumerate(refnas):
        if na.is_gap():
            continue
        bp = (bp + 1) % 3
        if bp == 1:
            bp1_indices.add(idx)

    if PositionalSeqStr.list_contains_any_gap(refnas):
        refnas = find_best_matches(refnas, seqnas, bp1_indices,
                                   scanstart=3, scanstep=3)
    if PositionalSeqStr.list_contains_any_gap(seqnas):
        seqnas = find_best_matches(seqnas, refnas, bp1_indices,
                                   scanstart=0, scanstep=1)
    return refnas, seqnas


NOGAP = 0b00
REFGAP = 0b01
SEQGAP = 0b10


def codon_pairs_group_key(cdpair):
    _, (refcd, seqcd) = cdpair
    if PositionalSeqStr.list_contains_any_gap(refcd):
        return REFGAP
    if PositionalSeqStr.list_contains_any_gap(seqcd):
        return SEQGAP
    return NOGAP


LEFT = 0b00
RIGHT = 0b01


def ensure_no_gap_ext_codons(ref_codons, seq_codons, side):
    if side == LEFT:
        ref_codons.reverse()
        seq_codons.reverse()
    idx = -1
    for idx, (refcd, seqcd) in enumerate(zip(ref_codons, seq_codons)):
        for na in chain(refcd, seqcd):
            if na.is_gap():
                break
        else:
            continue
        break
    else:
        # no gap in ext codons
        idx += 1
    ref_codons = ref_codons[:idx]
    seq_codons = seq_codons[:idx]
    if side == LEFT:
        ref_codons.reverse()
        seq_codons.reverse()
    return tuple(ref_codons), tuple(seq_codons), len(ref_codons)


def realign_gaps(refnas, seqnas, window_size):
    refcodons, seqcodons = group_by_codons(refnas, seqnas)
    refcodons = move_gap_to_codon_end(refcodons)
    seqcodons = move_gap_to_codon_end(seqcodons)

    # First gather gaps together according to window
    for pointer in range(0, len(refcodons) - window_size):
        slicekey = slice(pointer, pointer + window_size)
        win_refcodons = refcodons[slicekey]
        win_seqcodons = seqcodons[slicekey]
        win_refnas = list(chain(*win_refcodons))
        win_seqnas = list(chain(*win_seqcodons))
        win_refnas = window_gather_nearby_gaps(win_refnas)
        if PositionalSeqStr.list_contains_any_gap(win_refnas):
            # align gaps in seqnas with refnas
            # ref: AAA---CCCTTT     AAA---CCCTTT
            #                    =>
            # seq: AAACCC---TTT     AAA---CCCTTT
            gapidx = list_index(win_refnas, lambda na: na.is_gap())
            win_seqnas = move_gaps_to_index(win_seqnas, gapidx)
        else:
            win_seqnas = window_gather_nearby_gaps(win_seqnas)
        win_refnas, win_seqnas = clean_nalist_pairs(win_refnas, win_seqnas)
        (win_refcodons,
         win_seqcodons) = group_by_codons(win_refnas, win_seqnas)
        refcodons[slicekey] = win_refcodons
        seqcodons[slicekey] = win_seqcodons

    # Then adjust continuous refgap/seqgap placement
    gap_groups = groupby(
        enumerate(zip(refcodons,
                      seqcodons)),
        codon_pairs_group_key)
    for gap_type, codonpairs in gap_groups:
        if gap_type == NOGAP:
            continue
        refpos0, codonpairs = list(zip(*codonpairs))
        start, end = refpos0[0], refpos0[-1] + 1
        (refcds,
         seqcds) = list(zip(*codonpairs))

        # extend refcds/seqcds
        ext_refcds, ext_seqcds, offset = ensure_no_gap_ext_codons(
            refcodons[max(0, start - window_size):start],
            seqcodons[max(0, start - window_size):start],
            LEFT
        )
        refcds = ext_refcds + refcds
        seqcds = ext_seqcds + seqcds
        start -= offset

        ext_refcds, ext_seqcds, offset = ensure_no_gap_ext_codons(
            refcodons[end:end + window_size],
            seqcodons[end:end + window_size],
            RIGHT
        )
        refcds = refcds + ext_refcds
        seqcds = seqcds + ext_seqcds
        end += offset

        win_refnas = list(chain(*refcds))
        win_seqnas = list(chain(*seqcds))
        win_refnas, win_seqnas = \
            paired_find_best_matches(win_refnas, win_seqnas)
        (win_refcodons,
         win_seqcodons) = group_by_codons(win_refnas, win_seqnas)
        refcodons[start:end] = win_refcodons
        seqcodons[start:end] = win_seqcodons

    return list(chain(*refcodons)), list(chain(*seqcodons))


def move_gaps_to_index(nas, index):
    new_nas = [na for na in nas if not na.is_gap()]
    new_nas[index:index] = [na for na in nas if na.is_gap()]
    return new_nas


def window_gather_nearby_gaps(nas):
    if PositionalSeqStr.list_contains_any_gap(nas):
        gapidx = list_index(nas, lambda na: na.is_gap())
        return move_gaps_to_index(nas, gapidx)
    else:
        return nas


def gather_nearby_gaps(nas, window_size):
    na_winsize = window_size * 3
    for na_pointer in range(0, len(nas) - na_winsize):
        slicekey = slice(na_pointer, na_pointer + na_winsize)
        nas[slicekey] = window_gather_nearby_gaps(nas[slicekey])
    return nas


def clean_nalist_pairs(refnas, seqnas):
    """Remove matched gaps

    Example of matched gaps:

    TCAGCTGATGCA---CAA
    TCAGCTGATGCAC---AA
                 ^^
                 These two are "matched gaps" and can be
                 removed from pairwise alignment
    """
    return list(zip(*[
        (refna, seqna)
        for refna, seqna in zip(refnas, seqnas)
        if not refna.is_gap() or not seqna.is_gap()
    ]))


def codon_align(
    refseq,
    seq,
    window_size,
    refstart,
    refend,
    check_boundary=True
):

    reftext = refseq.seqtext
    if check_boundary and (reftext[0].is_gap() or reftext[-1].is_gap()):
        raise click.ClickException(
            'Unable to perform codon-alignment without the alignments '
            'being trimmed properly. Can be solved by pre-processing '
            'the alignments by command "trim-by-ref".'
        )

    idxstart, idxend = reftext.posrange2indexrange(refstart, refend)
    if idxstart == idxend:
        # nothing to be codon aligned
        return refseq, seq

    # step 1: apply reading frame
    reftext = reftext[idxstart:idxend]
    seqtext = seq.seqtext[idxstart:idxend]

    # step 2: remove matched gaps (due to MSA) from ref and seq
    refnas, seqnas = clean_nalist_pairs(reftext, seqtext)

    # step 3: gather and re-align nearby gaps located in same window
    refnas, seqnas = realign_gaps(refnas, seqnas, window_size)

    # last step: save "codon aligned" refseq and seq
    refseq = refseq.push_seqtext(
        refseq.seqtext[:idxstart] +
        PositionalSeqStr.join(refnas) +
        refseq.seqtext[idxend:],
        'codonalign({},{})'.format(refstart, refend), 0)
    seq = seq.push_seqtext(
        seq.seqtext[:idxstart] +
        PositionalSeqStr.join(seqnas) +
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
) -> Processor[Iterable[RefSeqPair]]:
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
        for refseq, seq in iterator:
            my_ref_end = ref_end
            if my_ref_end <= 0:
                my_ref_end = refseq.seqtext.max_pos
            if refseq.seqtype != 'NA':
                raise click.ClickException(
                    'Codon alignment only applies to nucleotide '
                    'sequences.')
            if len(seq.seqtext) == 0:
                # skip empty sequences
                yield refseq, seq
            else:
                yield codon_align(
                    refseq,
                    seq,
                    window_size,
                    ref_start, my_ref_end
                )

    return processor
