import re
from typing import TypeVar, List, Tuple

from ..models import Position

CIGAR_PATTERN = re.compile(r'(\d+)([MNDI])')

T = TypeVar('T', bound='CIGAR')


class CIGAR:

    ref_start: int
    seq_start: int
    cigar_string: str
    cigar_tuple: List[Tuple[int, str]]

    def __init__(
        self: T,
        ref_start: int,
        seq_start: int,
        cigar_string: str
    ) -> None:
        self.ref_start = ref_start
        self.seq_start = seq_start
        self.cigar_string = cigar_string
        self.cigar_tuple = [
            (int(num), op)
            for num, op in CIGAR_PATTERN.findall(cigar_string)
        ]

    def get_alignment(
        self: T,
        refseq: Position,
        seq: Position
    ) -> Tuple[Position, Position]:
        num: int
        op: str
        aligned_refseq: Position = refseq[self.ref_start:]
        aligned_seq: Position = seq[self.seq_start:]
        offset: int = 0
        for num, op in self.cigar_tuple:
            if op == 'M':
                offset += num
            elif op in ('D', 'N'):
                aligned_seq = (
                    aligned_seq[:offset] +
                    type(aligned_seq).init_gaps(gaplen=num) +
                    aligned_seq[offset:])
                offset += num
            elif op == 'I':
                aligned_refseq = (
                    aligned_refseq[:offset] +
                    type(aligned_refseq).init_gaps(gaplen=num) +
                    aligned_refseq[offset:])
                offset += num
        aligned_seq = aligned_seq[:offset]
        aligned_refseq = aligned_refseq[:offset]
        if len(aligned_refseq) != len(aligned_seq):
            raise ValueError(
                'Unmatched alignment length: {!r} and {!r}'
                .format(aligned_refseq, aligned_seq)
            )
        return aligned_refseq, aligned_seq

    def __repr__(self: T) -> str:
        return ('<CIGAR {!r} ref_start={!r} seq_start={!r}>'
                .format(self.cigar_string,
                        self.ref_start,
                        self.seq_start))
