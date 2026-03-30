from collections.abc import Generator, Iterable
from itertools import tee
from typing import TextIO

import click

from ..models import Position, RefSeqPair, Sequence
from . import fasta


def load(msafp: TextIO, reference: str, seqtype: type[Position]) -> Generator[RefSeqPair]:
    ref_finder: Iterable[Sequence]
    sequences: Iterable[Sequence] = fasta.load(msafp, seqtype)
    ref_finder, sequences = tee(sequences, 2)

    if reference:
        try:
            refseq = next(ref for ref in ref_finder if ref.header == reference)
        except StopIteration:
            raise click.ClickException(f'Unable to locate reference {reference!r} (--reference)') from None
    else:
        refseq = next(ref_finder)

    for seq in sequences:
        if seq != refseq:
            yield refseq, seq
