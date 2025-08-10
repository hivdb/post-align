import re
from typing import Any
from collections.abc import Iterable

from ..cli import cli
from ..models import Sequence, RefSeqPair

from ..processor import Processor, intermediate_processor

LEFT_TRIM_PATTERN: re.Pattern = re.compile(r'^[.-]+')
RIGHT_TRIM_PATTERN: re.Pattern = re.compile(r'[.-]+$')


def find_trim_slice(seq: Sequence) -> slice:
    left_trim: int | None = None
    seqtext_str: str = seq.seqtext_as_str
    match: re.Match | None = LEFT_TRIM_PATTERN.search(seqtext_str)
    if match:
        left_trim = match.span()[1]
    right_trim: int | None = None
    match = RIGHT_TRIM_PATTERN.search(seqtext_str)
    if match:
        right_trim = match.span()[0]
    return slice(left_trim, right_trim)


@cli.command('trim-by-ref')
def trim_by_ref() -> Processor[Iterable[RefSeqPair]]:
    """Trim all alignments by reference sequence"""

    @intermediate_processor('trim-by-ref')
    def processor(
        iterator: Iterable[RefSeqPair],
        *args: Any
    ) -> Iterable[RefSeqPair]:
        for refseq, seq in iterator:
            if seq.seqtext == '':
                # skip unaligned sequence
                yield refseq, seq
                continue
            sliceobj = find_trim_slice(refseq)
            refseq = refseq[sliceobj]
            seq = seq[sliceobj]
            yield refseq, seq

    return processor
