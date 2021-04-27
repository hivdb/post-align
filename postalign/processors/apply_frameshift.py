import re
import click
from functools import reduce
from operator import add
from more_itertools import chunked

from ..cli import cli


FS_PATTERN = re.compile(r'^(\d+)([+-]\d+)$')


def frameshift_callback(ctx, param, value):
    frameshift = {}
    for fs in value.split(','):
        match = FS_PATTERN.match(fs)
        if not match:
            raise click.BadArgumentUsage(
                'Invalid frameshift format: {!r}'.format(fs), ctx
            )
        pos, shift = match.groups()
        pos = int(pos, 10)
        shift = int(shift, 10)
        if not pos:
            raise click.BadArgumentUsage('POS cannot be zero')
        if not shift:
            raise click.BadArgumentUsage('N cannot be zero')
        if pos in frameshift:
            raise click.BadArgumentUsage(
                'POS {} is mentioned twice or more'
                .format(pos))
        frameshift[pos] = shift
    if not frameshift:
        raise click.BadArgumentUsage('argument FRAMESHIFT is required')
    return sorted(frameshift.items())


def apply_frameshift_to_single_seq(refseq, seq, frameshift):
    if len(seq) == 0:
        return refseq, seq

    breakpoints = [1]
    for pos, shift in frameshift:
        breakpoints.append(pos)
        breakpoints.append(pos + shift + 1)
    breakpoints.append(
        max(breakpoints[-1], refseq.seqtext.max_pos)
    )
    breakpoints = list(chunked(breakpoints, 2))
    adjusted_refseq = []
    adjusted_seq = []
    for pos_start, pos_end in breakpoints:
        idx_start, idx_end = refseq.seqtext.posrange2indexrange(
            pos_start, pos_end
        )
        adjusted_refseq.append(refseq[idx_start:idx_end])
        adjusted_seq.append(seq[idx_start:idx_end])
    adjusted_refseq = reduce(add, adjusted_refseq)
    adjusted_seq = reduce(add, adjusted_seq)

    return adjusted_refseq, adjusted_seq


@cli.command('apply-frameshift')
@click.argument('frameshift', callback=frameshift_callback)
def apply_frameshift(frameshift):
    """Apply known frameshifts to all alignments

    FRAMESHIFT (comma-delimited) Format: <POS1><+N|-N>,<POS2><+N|-N>...

    Where,
      POS: position in reference sequence, start from 1;
      +N: remove N notations started from given position;
      -N: repeat N notations started from given position.

    For example, in SARS-CoV-2 RdRP, between ORF1a/1b there is a ribosomal
    frameshift site (position 27) which pratically caused position 27 (a
    Cytosine, C) repeated once to produce a codon "CGG" as the amino acid
    postion R10.

      TCA GCT GAT GCA CAA TCG TTT TTA AAC .GG G...
                  1     1      2      2 2 2   3
           5      0     5      0      5 7 7   0
                                          ^
                                          repeat 27C

    To apply this known frameshift to all alignments, suppose the reference
    sequence started from the 'TCAGCT...' as above, a FRAMESHIFT argument
    should be formated as:

      27-1

    """

    def processor(iterator):

        for refseq, seq in iterator:
            yield apply_frameshift_to_single_seq(refseq, seq, frameshift)

    processor.command_name = 'apply-frameshift'
    processor.is_output_command = False
    return processor
