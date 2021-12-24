import click
from typing import TextIO, Generator, Iterable
from itertools import tee

from ..models import RefSeqPair, Sequence

from . import fasta


def load(
    msafp: TextIO,
    reference: str,
    seqtype: str
) -> Generator[RefSeqPair, None, None]:
    ref_finder: Iterable[Sequence]
    sequences: Iterable[Sequence] = fasta.load(msafp, seqtype)
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
