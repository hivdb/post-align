import re

import cython  # type: ignore

from ..models import NAPosition
from ..models.sequence import Position

CIGAR_PATTERN = re.compile(r'(\d+)([MNDI])')


@cython.ccall
@cython.returns(tuple)
def _get_alignment(
    cigar: 'CIGAR', refseq: list[NAPosition], seq: list[NAPosition], seqtype: type[Position]
) -> tuple[list[NAPosition], list[NAPosition]]:
    num: int
    op: str
    aligned_refseq: list[NAPosition] = list(refseq[cigar.ref_start :])
    aligned_seq: list[NAPosition] = list(seq[cigar.seq_start :])
    offset: int = 0
    for num, op in cigar.cigar_tuple:
        if op == 'M':
            offset += num
        elif op in ('D', 'N'):
            aligned_seq[offset:offset] = NAPosition.init_gaps(num)
            offset += num
        elif op == 'I':
            aligned_refseq[offset:offset] = NAPosition.init_gaps(num)
            offset += num
    aligned_seq = aligned_seq[:offset]
    aligned_refseq = aligned_refseq[:offset]
    if len(aligned_refseq) != len(aligned_seq):
        raise ValueError(
            f'Unmatched alignment length: {NAPosition.as_str(aligned_refseq)!r} and {NAPosition.as_str(aligned_seq)!r}'
        )

    return aligned_refseq, aligned_seq


@cython.cclass
class CIGAR:
    ref_start: int
    seq_start: int
    cigar_string: str
    cigar_tuple: list[tuple[int, str]]

    def __init__(self: 'CIGAR', ref_start: int, seq_start: int, cigar_string: str) -> None:
        self.ref_start = ref_start
        self.seq_start = seq_start
        self.cigar_string = cigar_string
        self.cigar_tuple = [(int(num), op) for num, op in CIGAR_PATTERN.findall(cigar_string)]

    def get_cigar_string(self: 'CIGAR') -> str:
        return self.cigar_string

    def shrink_by_ref(self: 'CIGAR', keep_size: int) -> 'CIGAR':
        new_string: str
        new_tuple: list[tuple[int, str]] = []
        remain_size: int = keep_size
        for num, op in self.cigar_tuple:
            if op == 'I':
                new_tuple.append((num, op))
            elif num < remain_size:
                new_tuple.append((num, op))
                remain_size -= num
            else:
                new_tuple.append((remain_size, op))
                break
        new_string = ''.join([f'{num}{op}' for num, op in new_tuple])
        return CIGAR(self.ref_start, self.seq_start, new_string)

    def get_alignment(
        self: 'CIGAR',
        refseq: list[NAPosition],
        seq: list[NAPosition],
        seqtype: type[Position],
    ) -> tuple[list[NAPosition], list[NAPosition]]:
        alignment: tuple[list[NAPosition], list[NAPosition]] = _get_alignment(self, refseq, seq, seqtype)
        return alignment

    def __repr__(self: 'CIGAR') -> str:
        return f'<CIGAR {self.cigar_string!r} ref_start={self.ref_start!r} seq_start={self.seq_start!r}>'
