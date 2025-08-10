import cython  # type: ignore
from typing import Any
from collections.abc import ByteString
from itertools import zip_longest
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
    nas: list['NAPosition'],
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
    nas: list['NAPosition'],
    pos_start: int,
    pos_end: int,
    include_boundary_gaps: bool = False
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
        payload: Any = None
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
        return (
            'NAPosition({!r}, {!r}, {!r}, {!r})'
            .format(self.notation, self.pos, self.flag, self.payload)
        )

    def __copy__(self: 'NAPosition') -> 'NAPosition':
        return NAPosition(self.notation, self.pos, self.flag, self.payload)

    @classmethod
    def init_gaps(
        cls: type['NAPosition'],
        gaplen: int
    ) -> list['NAPosition']:
        return [
            cls(GAP_CHAR, -1, PositionFlag.NONE)
            for _ in range(gaplen)
        ]

    @classmethod
    def init_from_bytes(
        cls: type['NAPosition'],
        seq_text: ByteString,
        seq_payload: list | None = None
    ) -> list['NAPosition']:
        """Generate positions from raw sequence bytes.

        :param seq_text: Raw nucleotide sequence data.
        :param seq_payload: Optional per-position payload.
        :returns: List of initialized positions.
        """
        seq_text = bytes(seq_text).upper()
        if seq_payload is None:
            seq_payload = []
        return [
            cls(na, pos, PositionFlag.NONE, payload)
            for na, pos, payload in zip_longest(
                seq_text,
                enumerate_seq_pos(seq_text),
                seq_payload
            )
        ]

    @staticmethod
    def min_pos(nas: list['NAPosition']) -> int:
        """Return the smallest positive position.

        :param nas: Positions to inspect.
        :returns: Minimum non-zero position or ``-1`` if none found.
        """
        for na in nas:
            if na.pos > 0:
                return na.pos
        return -1

    @staticmethod
    def max_pos(nas: list['NAPosition']) -> int:
        """Return the largest positive position.

        :param nas: Positions to inspect.
        :returns: Maximum non-zero position or ``-1`` if none found.
        """
        for na in reversed(nas):
            if na.pos > 0:
                return na.pos
        return -1

    @staticmethod
    def min_nongap_index(
        nas: list['NAPosition'],
        start: int = -1,
        stop: int = -1
    ) -> int:
        """Index of first non-gap position.

        :param nas: Positions to inspect.
        :param start: Optional starting index.
        :param stop: Optional stopping index.
        :returns: Index of first non-gap or ``-1`` if none found.
        """
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
        nas: list['NAPosition'],
        start: int = -1,
        stop: int = -1
    ) -> int:
        """Index of last non-gap position.

        :param nas: Positions to inspect.
        :param start: Optional starting index.
        :param stop: Optional stopping index.
        :returns: Index of last non-gap or ``-1`` if none found.
        """
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
        """Apply a flag to all positions.

        :param nas: Positions to modify.
        :param flag: Flag to set.
        """
        for na in nas:
            na.flag |= flag

    @staticmethod
    def any_has_flag(nas: list['NAPosition'], flag: PositionFlag) -> bool:
        """Check if any position has a specific flag.

        :param nas: Positions to inspect.
        :param flag: Flag to test.
        :returns: ``True`` if any position has the flag.
        """
        return any([
            flag & na.flag for na in nas
        ])

    @staticmethod
    def all_have_flag(nas: list['NAPosition'], flag: PositionFlag) -> bool:
        """Check if all positions have a specific flag.

        :param nas: Positions to inspect.
        :param flag: Flag to test.
        :returns: ``True`` if all positions have the flag.
        """
        return all([
            flag & na.flag for na in nas
        ])

    posrange2indexrange = _posrange2indexrange

    @staticmethod
    def count_gaps(nas: list['NAPosition']) -> int:
        """Count gap positions.

        :param nas: Positions to inspect.
        :returns: Number of gaps.
        """
        return sum([na.is_gap for na in nas])

    @staticmethod
    def count_nongaps(nas: list['NAPosition']) -> int:
        """Count non-gap positions.

        :param nas: Positions to inspect.
        :returns: Number of non-gaps.
        """
        return sum([not na.is_gap for na in nas])

    @staticmethod
    def remove_gaps(
        nas: list['NAPosition']
    ) -> list['NAPosition']:
        """Filter out gap positions.

        :param nas: Positions to inspect.
        :returns: List without gap positions.
        """
        return [na for na in nas if not na.is_gap]

    @staticmethod
    def as_bytes(
        nas: list['NAPosition']
    ) -> bytes:
        """Serialize positions as raw bytes.

        :param nas: Positions to convert.
        :returns: Sequence encoded as bytes.
        """
        return bytes([na.notation for na in nas])

    @classmethod
    def as_str(
        cls: type['NAPosition'],
        nas: list['NAPosition']
    ) -> str:
        return str(cls.as_bytes(nas), 'ASCII')

    @staticmethod
    def any_has_gap(nas: list['NAPosition']) -> bool:
        """Check if any position is a gap.

        :param nas: Positions to inspect.
        :returns: ``True`` if any gap is present.
        """
        return any([na.is_gap for na in nas])

    @staticmethod
    def all_have_gap(nas: list['NAPosition']) -> bool:
        """Check if all positions are gaps.

        :param nas: Positions to inspect.
        :returns: ``True`` if all are gaps.
        """
        return all([na.is_gap for na in nas])
