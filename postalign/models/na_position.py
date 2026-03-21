from collections.abc import ByteString
from itertools import zip_longest
from typing import Any

import cython  # type: ignore

from .position_flag import PositionFlag

GAP_CHAR: int = ord(b'-')
GAP_CHARS: tuple[int, ...] = tuple(b'-.')

FIRST: cython.int = 0
LAST: cython.int = 1


@cython.ccall
@cython.returns(list)
def enumerate_seq_pos(seq_text: bytes) -> list[int]:
    na: int
    offset: int = 1
    seq_pos: list[int] = []
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
    nas: list['NAPosition'], pos: cython.int, first: cython.int
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
    nas: list['NAPosition'],
    pos_start: int,
    pos_end: int,
    include_boundary_gaps: bool = False,
) -> tuple[int, int]:
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
        pos_start = min_pos if min_pos > pos_start else pos_start
        pos_end = max_pos if max_pos < pos_end else pos_end
        for pos in range(pos_start, pos_end + 1):
            idx_start = _pos2index(nas, pos, FIRST)
            if idx_start > -1:
                break
        for pos in range(pos_end, pos_start - 1, -1):
            idx_end = _pos2index(nas, pos, LAST) + 1
            if idx_end > -1:
                break
    if include_boundary_gaps:
        while idx_start > 0 and nas[idx_start - 1].is_gap:
            idx_start -= 1
        naslen: int = len(nas)
        while idx_end < naslen and nas[idx_end].is_gap:
            idx_end += 1
    return idx_start, idx_end


@cython.cclass
class NAPosition:
    notation: int = cython.declare(cython.int, visibility='public')
    pos: int = cython.declare(cython.int, visibility='public')
    flag: PositionFlag = cython.declare(cython.int, visibility='public')
    is_gap: bool = cython.declare(cython.bint, visibility='public')

    # `payload` can carry any kind of object. It is useful in circumstances
    # such as codon-aligning multi-alignment sequences
    payload: Any = cython.declare(object, visibility='public')

    def __init__(
        self: 'NAPosition',
        notation: int,
        pos: int,
        flag: PositionFlag,
        payload: Any = None,
    ) -> None:
        self.notation = notation
        self.pos = pos
        self.flag = flag
        self.is_gap = notation in GAP_CHARS
        self.payload = payload

    def __str__(self: 'NAPosition') -> str:
        return chr(self.notation)

    def __bytes__(self: 'NAPosition') -> bytes:
        return bytes([self.notation])

    def __repr__(self: 'NAPosition') -> str:
        return f'NAPosition({self.notation!r}, {self.pos!r}, {self.flag!r}, {self.payload!r})'

    def __copy__(self: 'NAPosition') -> 'NAPosition':
        return NAPosition(self.notation, self.pos, self.flag, self.payload)

    @classmethod
    def init_gaps(cls: type['NAPosition'], gaplen: int) -> list['NAPosition']:
        return [cls(GAP_CHAR, -1, PositionFlag.NONE) for _ in range(gaplen)]

    @classmethod
    def init_from_bytes(
        cls: type['NAPosition'],
        seq_text: ByteString,
        seq_payload: list | None = None,
    ) -> list['NAPosition']:
        na: int  # noqa: F842
        pos: int  # noqa: F842
        seq_text = bytes(seq_text).upper()
        if seq_payload is None:
            seq_payload = []
        return [
            cls(na, pos, PositionFlag.NONE, payload)
            for na, pos, payload in zip_longest(
                seq_text, enumerate_seq_pos(seq_text), seq_payload
            )
        ]

    @staticmethod
    def min_pos(nas: list['NAPosition']) -> int:
        na: NAPosition
        for na in nas:
            if na.pos > 0:
                return na.pos
        return -1

    @staticmethod
    def max_pos(nas: list['NAPosition']) -> int:
        na: NAPosition
        for na in reversed(nas):
            if na.pos > 0:
                return na.pos
        return -1

    @staticmethod
    def min_nongap_index(
        nas: list['NAPosition'], start: int = -1, stop: int = -1
    ) -> int:
        na: NAPosition
        for idx, na in enumerate(nas):
            if start > -1 and idx < start:
                continue
            if stop > -1 and idx >= stop:
                break
            if na.pos > 0:
                return idx
        return -1

    @staticmethod
    def max_nongap_index(
        nas: list['NAPosition'], start: int = -1, stop: int = -1
    ) -> int:
        na: NAPosition
        length = len(nas)
        for revidx, na in enumerate(reversed(nas)):
            idx = length - 1 - revidx
            if start > -1 and idx < start:
                break
            if stop > -1 and idx >= stop:
                continue
            if na.pos > 0:
                return idx
        return -1

    @staticmethod
    def set_flag(nas: list['NAPosition'], flag: PositionFlag) -> None:
        na: NAPosition
        for na in nas:
            na.flag |= flag

    @staticmethod
    def any_has_flag(nas: list['NAPosition'], flag: PositionFlag) -> bool:
        na: NAPosition  # noqa: F842
        return any(flag & na.flag for na in nas)

    @staticmethod
    def all_have_flag(nas: list['NAPosition'], flag: PositionFlag) -> bool:
        na: NAPosition  # noqa: F842
        return all(flag & na.flag for na in nas)

    posrange2indexrange = _posrange2indexrange

    @staticmethod
    def count_gaps(nas: list['NAPosition']) -> int:
        na: NAPosition  # noqa: F842
        return sum(na.is_gap for na in nas)

    @staticmethod
    def count_nongaps(nas: list['NAPosition']) -> int:
        na: NAPosition  # noqa: F842
        return sum(not na.is_gap for na in nas)

    @staticmethod
    def remove_gaps(nas: list['NAPosition']) -> list['NAPosition']:
        na: NAPosition  # noqa: F842
        return [na for na in nas if not na.is_gap]

    @staticmethod
    def as_bytes(nas: list['NAPosition']) -> bytes:
        na: NAPosition  # noqa: F842
        return bytes([na.notation for na in nas])

    @classmethod
    def as_str(cls: type['NAPosition'], nas: list['NAPosition']) -> str:
        return str(cls.as_bytes(nas), 'ASCII')

    @staticmethod
    def any_has_gap(nas: list['NAPosition']) -> bool:
        na: NAPosition  # noqa: F842
        return any(na.is_gap for na in nas)

    @staticmethod
    def all_have_gap(nas: list['NAPosition']) -> bool:
        na: NAPosition  # noqa: F842
        return all(na.is_gap for na in nas)
