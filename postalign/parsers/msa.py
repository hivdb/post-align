import click
from itertools import tee

from . import fasta


def load(msafp, reference, seqtype):
    sequences = fasta.load(msafp, seqtype)
    ref_finder, sequences = tee(sequences, 2)

    if reference:
        try:
            refseq = next(
                ref
                for ref in ref_finder
                if ref.header == reference
            )
        except StopIteration:
            raise click.ClickException(
                'Unable to locate reference {!r} (--reference)'
                .format(reference)
            )
    else:
        refseq = next(ref_finder)

    for seq in sequences:
        if seq != refseq:
            yield refseq, seq
