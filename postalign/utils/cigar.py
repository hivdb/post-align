import re
import cython  # type: ignore
from typing import TypeVar, List, Tuple, Type

from ..models import Position

CIGAR_PATTERN = re.compile(r'(\d+)([MNDI])')

T = TypeVar('T', bound='CIGAR')


@cython.ccall
@cython.returns(tuple)
def _get_alignment(
    cigar: 'CIGAR',
    refseq: List[Position],
    seq: List[Position],
    seqtype: Type[Position]
) -> Tuple[List[Position], List[Position]]:
    num: int
    op: str
    aligned_refseq: List[Position] = refseq[cigar.ref_start:]
    aligned_seq: List[Position] = seq[cigar.seq_start:]
    offset: int = 0
    for num, op in cigar.cigar_tuple:
        if op == 'M':
            offset += num
        elif op in ('D', 'N'):
            aligned_seq[offset:offset] = seqtype.init_gaps(num)
            offset += num
        elif op == 'I':
            aligned_refseq[offset:offset] = seqtype.init_gaps(num)
            offset += num
    aligned_seq = aligned_seq[:offset]
    aligned_refseq = aligned_refseq[:offset]
    if len(aligned_refseq) != len(aligned_seq):
        raise ValueError(
            'Unmatched alignment length: {!r} and {!r}'
            .format(
                seqtype.as_str(aligned_refseq),
                seqtype.as_str(aligned_seq)
            )
        )

    return aligned_refseq, aligned_seq


@cython.cclass
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

    def get_cigar_string(self: T) -> str:
        return self.cigar_string

    def shrink_by_ref(self: T, keep_size: int):
        new_string: str
        new_tuple: List[Tuple[int, str]] = []
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
        new_string = ''.join(['{}{}'.format(num, op) for num, op in new_tuple])
        return CIGAR(
            self.ref_start,
            self.seq_start,
            new_string
        )

    def get_alignment(
        self: T,
        refseq: List[Position],
        seq: List[Position],
        seqtype: Type[Position]
    ) -> Tuple[List[Position], List[Position]]:
        alignment: Tuple[
            List[Position],
            List[Position]
        ] = _get_alignment(self, refseq, seq, seqtype)
        return alignment

    def __repr__(self: T) -> str:
        return ('<CIGAR {!r} ref_start={!r} seq_start={!r}>'
                .format(self.cigar_string,
                        self.ref_start,
                        self.seq_start))
