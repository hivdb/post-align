import click
from collections import namedtuple

from ..cli import cli


GAP_CHARS = '.-'
NAPair = namedtuple(
    'NAPair', ['refna', 'seqna', 'refgap', 'seqgap', 'refbp']
)


def list_index(listdata, cond, start=0):
    for idx, n in enumerate(listdata[start:]):
        if cond(n):
            return start + idx, n
    raise ValueError(
        'given condition {!r} is not found in list'
        .format(cond)
    )


def list_rindex(listdata, cond, start=None):
    if start is None:
        start = len(listdata)
    for idx, n in enumerate(listdata[:start][::-1]):
        if cond(n):
            return start - idx - 1, n
    raise ValueError(
        'given condition {!r} is not found in list'
        .format(cond)
    )


def combine_nearby_gaps(napairs, naattr, gapattr):

    def getgap(np):
        return getattr(np, gapattr)

    def getna(np):
        return getattr(np, naattr)

    gapidx = list_index(napairs, getgap)
    nas = [
        getna(np) for np in napairs if getgap(np)
    ]
    nas[gapidx:gapidx] = [
        getna(np) for np in napairs if not getgap(np)
    ]
    return [
        NAPair(**{
            **np,
            naattr: na
        })
        for na, np in zip(nas, napairs)
    ]


def codon_align(refseq, seq, reading_frame, window_size):
    reftext = refseq.seqtext[reading_frame - 1:]
    seqtext = seq.seqtext[reading_frame - 1:]
    napairs = []
    bp = 2
    for refna, seqna in zip(reftext, seqtext):
        ref_is_gap = refna in GAP_CHARS
        seq_is_gap = seqna in GAP_CHARS
        if ref_is_gap and seq_is_gap:
            continue
        if not ref_is_gap:
            bp = (bp + 1) % 3
        napairs.append(NAPair(
            refna, seqna, int(ref_is_gap), int(seq_is_gap), bp
        ))

    na_winsize = window_size * 3
    for na_pointer in range(0, len(napairs) - na_winsize):
        slicekey = slice(na_pointer, na_pointer + na_winsize)
        win_napairs = napairs[slicekey]
        len_refgaps = sum(np.refgap for np in win_napairs)
        len_seqgaps = sum(np.seqgap for np in win_napairs)
        if len_refgaps == len_seqgaps == 0:
            # no gap exists in current window, move on
            continue
        if len_refgaps > 1:
            win_napairs = combine_nearby_gaps(win_napairs, 'refna', 'refgap')
        if len_seqgaps > 1:
            win_napairs = combine_nearby_gaps(win_napairs, 'seqna', 'seqgap')
        napairs[slicekey] = win_napairs

    reftext = ''.join(np.refna for np in napairs)
    seqtext = ''.join(np.seqna for np in napairs)
    refseq = refseq.modify_seqtext(reftext, 'codonalign()')
    seq = seq.modify_seqtext(seqtext, 'codonalign()')
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
            yield codon_align(refseq, seq, reading_frame, window_size)

    processor.command_name = 'codon-alignment'
    processor.is_output_command = False
    return processor
