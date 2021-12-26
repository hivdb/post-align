import re
from typing import TypeVar, List, Tuple, Type

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
        refseq: List[Position],
        seq: List[Position],
        seqtype: Type[Position]
    ) -> Tuple[List[Position], List[Position]]:
        num: int
        op: str
        aligned_refseq: List[Position] = refseq[self.ref_start:]
        aligned_seq: List[Position] = seq[self.seq_start:]
        offset: int = 0
        for num, op in self.cigar_tuple:
            if op == 'M':
                offset += num
            elif op in ('D', 'N'):
                aligned_seq = (
                    aligned_seq[:offset] +
                    seqtype.init_gaps(num) +
                    aligned_seq[offset:])
                offset += num
            elif op == 'I':
                aligned_refseq = (
                    aligned_refseq[:offset] +
                    seqtype.init_gaps(num) +
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
