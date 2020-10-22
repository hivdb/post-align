import re
import click
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
        left_start = left_end = adjusted_refseq = None
        for refseq, seq in iterator:
            if len(seq) == 0:
                yield refseq, seq
                continue
            if adjusted_refseq is None:
                reflen = len(refseq)
                breakpoints = [0]
                for pos, shift in frameshift:
                    if pos > reflen:
                        # frameshift is longer than reference
                        break
                    breakpoints.append(pos)
                    breakpoints.append(pos + shift)
                breakpoints.append(reflen)
                [(left_start, left_end),
                 *slicetuples] = list(chunked(breakpoints, 2))
                adjusted_refseq = refseq[left_start:left_end]
                for start, end in slicetuples:
                    adjusted_refseq += refseq[start:end]
            adjusted_seq = seq[left_start:left_end]
            for start, end in slicetuples:
                adjusted_seq += seq[start:end]
            yield (adjusted_refseq, adjusted_seq)

    processor.command_name = 'apply-frameshift'
    processor.is_output_command = False
    return processor
