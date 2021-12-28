import cython  # type: ignore
from typing import (
    List,
    Tuple,
    ByteString,
    Type
)
from enum import IntFlag

GAP_CHAR: int = ord(b'-')
GAP_CHARS: Tuple[int, ...] = tuple(b'-.')

FIRST: cython.int = 0
LAST: cython.int = 1


class NAFlag(IntFlag):
    TRIM_BY_SEQ = 1
    NONE = 0


@cython.ccall
@cython.returns(list)
def enumerate_seq_pos(seq_text: bytes) -> List[int]:
    na: int
    offset: int = 1
    seq_pos: List[int] = []
    for na in seq_text:
        if na in GAP_CHARS:
            seq_pos.append(-1)
        else:
            seq_pos.append(offset)
            offset += 1
    return seq_pos


@cython.cfunc
@cython.inline
def _pos2index(
    nas: List['NAPosition'],
    pos: cython.int,
    first: cython.int
) -> cython.int:
    idx: cython.int
    len_nas: cython.int = len(nas)
    if first == FIRST:
        idx = 0
        while idx < len_nas:
            if nas[idx].pos == pos:
                return idx
            idx += 1
    elif first == LAST:
        idx = len_nas
        while idx > -1:
            idx -= 1
            if nas[idx].pos == pos:
                return idx
    return -1


@cython.ccall
@cython.returns(tuple)
def _posrange2indexrange(
    nas: List['NAPosition'],
    pos_start: int,
    pos_end: int
) -> Tuple[int, int]:
    idx_start: int
    idx_end: int
    max_pos: int = NAPosition.max_pos(nas)
    min_pos: int = NAPosition.min_pos(nas)
    if max_pos < 0 or min_pos < 0:
        return 0, 0
    if pos_start > max_pos:
        idx_start = idx_end = _pos2index(nas, max_pos, LAST)
    elif pos_end < min_pos:
        idx_start = idx_end = _pos2index(nas, min_pos, FIRST)
    else:
        pos_start = max(min_pos, pos_start)
        pos_end = min(max_pos, pos_end)
        for pos in range(pos_start, pos_end + 1):
            idx_start = _pos2index(nas, pos, FIRST)
            if idx_start > -1:
                break
        for pos in range(pos_end, pos_start - 1, -1):
            idx_end = _pos2index(nas, pos, LAST) + 1
            if idx_end > -1:
                break
    return idx_start, idx_end


@cython.cclass
class NAPosition:

    notation: int = cython.declare(cython.int, visibility='public')
    pos: int = cython.declare(cython.int, visibility='public')
    flag: NAFlag = cython.declare(cython.int, visibility='public')
    is_gap: bool = cython.declare(cython.bint, visibility='public')

    def __init__(
        self: 'NAPosition',
        notation: int,
        pos: int,
        flag: NAFlag
    ) -> None:
        self.notation = notation
        self.pos = pos
        self.flag = flag
        self.is_gap = notation in GAP_CHARS

    def __str__(self: 'NAPosition') -> str:
        return chr(self.notation)

    def __bytes__(self: 'NAPosition') -> bytes:
        return bytes([self.notation])

    def __repr__(self: 'NAPosition') -> str:
        return (
            'NAPosition({!r}, {!r}, {!r})'
            .format(self.notation, self.pos, self.flag)
        )

    @classmethod
    def init_gaps(
        cls: Type['NAPosition'],
        gaplen: int
    ) -> List['NAPosition']:
        return [
            cls(GAP_CHAR, -1, NAFlag.NONE)
            for _ in range(gaplen)
        ]

    @classmethod
    def init_from_triplets(
        cls: Type['NAPosition'],
        triplets: List[Tuple[int, int, NAFlag]]
    ) -> List['NAPosition']:
        return [cls(*triplet) for triplet in triplets]

    @classmethod
    def init_from_bytes(
        cls: Type['NAPosition'],
        seq_text: ByteString
    ) -> List['NAPosition']:
        na: int
        pos: int
        seq_text = bytes(seq_text).upper()
        return [
            cls(na, pos, NAFlag.NONE)
            for na, pos in zip(seq_text, enumerate_seq_pos(seq_text))
        ]

    @staticmethod
    def min_pos(nas: List['NAPosition']) -> int:
        na: NAPosition
        for na in nas:
            if na.pos > 0:
                return na.pos
        return -1

    @staticmethod
    def max_pos(nas: List['NAPosition']) -> int:
        na: NAPosition
        for na in reversed(nas):
            if na.pos > 0:
                return na.pos
        return -1

    @staticmethod
    def set_flag(nas: List['NAPosition'], flag: NAFlag) -> None:
        na: NAPosition
        for na in nas:
            na.flag |= flag

    @staticmethod
    def any_has_flag(nas: List['NAPosition'], flag: NAFlag) -> bool:
        na: NAPosition
        return any([
            flag & na.flag for na in nas
        ])

    @staticmethod
    def all_have_flag(nas: List['NAPosition'], flag: NAFlag) -> bool:
        na: NAPosition
        return all([
            flag & na.flag for na in nas
        ])

    posrange2indexrange = _posrange2indexrange

    @staticmethod
    def count_gaps(nas: List['NAPosition']) -> int:
        na: NAPosition
        return sum([na.is_gap for na in nas])

    @staticmethod
    def count_nongaps(nas: List['NAPosition']) -> int:
        na: NAPosition
        return sum([not na.is_gap for na in nas])

    @staticmethod
    def remove_gaps(
        nas: List['NAPosition']
    ) -> List['NAPosition']:
        na: NAPosition
        return [na for na in nas if not na.is_gap]

    @staticmethod
    def as_bytes(
        nas: List['NAPosition']
    ) -> bytes:
        na: NAPosition
        return bytes([na.notation for na in nas])

    @classmethod
    def as_str(
        cls: Type['NAPosition'],
        nas: List['NAPosition']
    ) -> str:
        return str(cls.as_bytes(nas), 'ASCII')

    @staticmethod
    def any_has_gap(nas: List['NAPosition']) -> bool:
        na: NAPosition
        return any([na.is_gap for na in nas])

    @staticmethod
    def all_have_gap(nas: List['NAPosition']) -> bool:
        na: NAPosition
        return all([na.is_gap for na in nas])
