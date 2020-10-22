import click
from itertools import chain, groupby

from ..cli import cli
from ..utils import group_by_codons


def list_index(listdata, cond, start=0):
    for idx, n in enumerate(listdata[start:]):
        if cond(n):
            return start + idx
    raise ValueError(
        'given condition {!r} is not found in list'
        .format(cond)
    )


def list_rindex(listdata, cond, start=None):
    if start is None:
        start = len(listdata)
    for idx, n in enumerate(listdata[:start][::-1]):
        if cond(n):
            return start - idx - 1
    raise ValueError(
        'given condition {!r} is not found in list'
        .format(cond)
    )


def move_gap_to_codon_end(codons):
    new_codons = []
    for codon in codons:
        new_codons.append(
            [na for na in codon if not na.is_gap()] +
            [na for na in codon if na.is_gap()]
        )
    return new_codons


def calc_match_score(mynas, othernas):
    # TODO: use blosum62 to calculate
    return sum(str(myna) == str(otherna)
               for myna, otherna in zip(mynas, othernas))


def find_best_matches(mynas, othernas, bp1_indices, scanstart=0, scanstep=1):
    orig_gapidx = list_index(mynas, lambda na: na.is_gap())
    mygap = [na for na in mynas if na.is_gap()]
    mynas = [na for na in mynas if not na.is_gap()]
    max_score = (-1, -1)
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
        if score > max_score:
            max_score = score
            best_mynas = test_mynas
    if best_mynas is None:
        best_mynas = mygap
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

    if any(na.is_gap() for na in refnas):
        refnas = find_best_matches(refnas, seqnas, bp1_indices,
                                   scanstart=3, scanstep=3)
    elif any(na.is_gap() for na in seqnas):
        seqnas = find_best_matches(seqnas, refnas, bp1_indices,
                                   scanstart=0, scanstep=1)
    return refnas, seqnas


NOGAP = 0b00
REFGAP = 0b01
SEQGAP = 0b10


def codon_pairs_group_key(cdpair):
    _, (refcd, seqcd) = cdpair
    if any(na.is_gap() for na in refcd):
        return REFGAP
    if any(na.is_gap() for na in seqcd):
        return SEQGAP
    return NOGAP


LEFT = 0b00
RIGHT = 0b01


def ensure_no_gap_ext_codons(ref_codons, seq_codons, side):
    if side == LEFT:
        ref_codons.reverse()
        seq_codons.reverse()
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
        if any(na.is_gap() for na in win_refnas):
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
    if any(na.is_gap() for na in nas):
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
    new_refnas = []
    new_seqnas = []
    for refna, seqna in zip(refnas, seqnas):
        if refna.is_gap() and seqna.is_gap():
            continue
        new_refnas.append(refna)
        new_seqnas.append(seqna)
    return new_refnas, new_seqnas


def codon_align(refseq, seq, reading_frame, window_size):

    reftext = refseq.seqtext
    if reftext[0].is_gap() or reftext[-1].is_gap():
        raise click.ClickException(
            'Unable to perform codon-alignment without the alignments '
            'being trimmed properly. Can be solved by pre-processing '
            'the alignments by command "trim-by-ref".'
        )

    # step 1: apply reading frame
    refstart = seqstart = reading_frame - 1
    reftext = reftext[refstart:]
    seqtext = seq.seqtext[seqstart:]

    # step 2: remove matched gaps (due to MSA) from ref and seq
    refnas, seqnas = clean_nalist_pairs(reftext, seqtext)

    # step 3: gather and re-align nearby gaps located in same window
    refnas, seqnas = realign_gaps(refnas, seqnas, window_size)

    # last step: save "codon aligned" refseq and seq
    refseq = refseq.push_seqtext(refnas, 'codonalign()', refstart)
    seq = seq.push_seqtext(seqnas, 'codonalign()', seqstart)
    return refseq, seq


@cli.command('codon-alignment')
@click.option(
    '--reading-frame',
    type=click.Choice([1, 2, 3]),
    default=1,
    help=(
        'Reading frame where codon-alignment started'
    ))
@click.option(
    '--window-size',
    type=int,
    default=5,
    help=(
        'Window size as # of codons for frameshift compensation'
    ))
def codon_alignment(reading_frame, window_size):
    """Codon alignment

    A blackbox re-implementation of the "codon-align" tool
    created by LANL HIV Sequence Database.
    """

    def processor(iterator):
        for refseq, seq in iterator:
            if refseq.seqtype != 'NA':
                raise click.ClickException(
                    'Codon alignment only applies to nucleotide '
                    'sequences.')
            if seq.seqtext == '':
                # skip empty sequences
                yield refseq, seq
            else:
                yield codon_align(refseq, seq, reading_frame, window_size)

    processor.command_name = 'codon-alignment'
    processor.is_output_command = False
    return processor
