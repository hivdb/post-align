"""MSA parser utilities."""

import typer
from typing import Type, TextIO, Generator, Iterable
from itertools import tee

from ..models import RefSeqPair, Sequence, Position

from . import fasta


def load(
    msafp: TextIO,
    reference: str,
    seqtype: Type[Position]
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
            raise typer.BadParameter(
                'Unable to locate reference {!r} (--reference)'.format(
                    reference
                )
            )
    else:
        refseq = next(ref_finder)

    for seq in sequences:
        if seq != refseq:
            yield refseq, seq
