import cython  # type: ignore
from typing import (
    List,
    Set,
    Optional,
    Tuple,
    Union,
    Iterable,
    Generator,
    ByteString,
    Type,
    TypeVar
)

GAP_CHARS = b'.-'


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


@cython.cclass
class NAPosition:

    @classmethod
    def init_empty(cls: Type['NAPosition']) -> 'NAPosition':
        return cls(b'', [], [])

    @classmethod
    def init_gaps(cls: Type['NAPosition'], *, gaplen: int) -> 'NAPosition':
        seq_text: bytes = b'-' * gaplen
        seq_pos: List[int] = [-1] * gaplen
        seq_flag: List[Set[str]] = [set() for _ in seq_pos]
        return cls(seq_text, seq_pos, seq_flag)

    @classmethod
    def init_from_triplets(
        cls: Type['NAPosition'],
        triplets: List[Tuple[int, int, Set[str]]]
    ) -> 'NAPosition':
        seq_text: bytearray = bytearray()
        seq_pos: List[int] = []
        seq_flag: List[Set[str]] = []
        for na, pos, flag in triplets:
            seq_text.append(na)
            seq_pos.append(pos)
            seq_flag.append(flag)
        return cls(bytes(seq_text), seq_pos, seq_flag)

    @classmethod
    def init_from_nastring(
        cls: Type['NAPosition'],
        seq_text: ByteString
    ) -> 'NAPosition':
        seq_text = bytes(seq_text).upper()
        seq_pos: List[int] = enumerate_seq_pos(seq_text)
        seq_flag: List[Set[str]] = [set() for _ in seq_pos]
        return cls(seq_text, seq_pos, seq_flag)

    _seq_text: bytes = cython.declare(bytes, visibility='public')
    _seq_pos: List[int] = cython.declare(list, visibility='public')
    _seq_flag: List[Set[str]] = cython.declare(list, visibility='public')
    is_single_gap: bool = cython.declare(cython.bint, visibility='public')
    _min_pos: int
    _max_pos: int
    _is_gap: Optional[bool]

    def __init__(
        self: 'NAPosition',
        seq_text: bytes,
        seq_pos: List[int],
        seq_flag: List[Set[str]]
    ) -> None:
        self._seq_text = seq_text
        self._seq_pos = seq_pos
        self._seq_flag = seq_flag
        self._min_pos = -1
        self._max_pos = -1
        if len(seq_text) == 1:
            self._is_gap = seq_text in GAP_CHARS
            self.is_single_gap = self._is_gap
        else:
            self._is_gap = None
            self.is_single_gap = False

    def min_pos(self: 'NAPosition') -> int:
        if self._min_pos < 0:
            self._min_pos = next(
                (pos for pos in self._seq_pos if pos > 0), -1
            )
        return self._min_pos

    def max_pos(self: 'NAPosition') -> int:
        if self._max_pos < 0:
            self._max_pos = next(
                (pos for pos in reversed(self._seq_pos) if pos > 0), -1
            )
        return self._max_pos

    def empty(self: 'NAPosition') -> bool:
        return self.min_pos() < 0

    def set_flag(self: 'NAPosition', flag: str) -> None:
        one: Set[str]
        for one in self._seq_flag:
            one.add(flag)

    def any_has_flag(self: 'NAPosition', flag: str) -> bool:
        one: Set[str]
        for one in self._seq_flag:
            if flag in one:
                return True
        return False

    def all_have_flag(self: 'NAPosition', flag: str) -> bool:
        one: Set[str]
        return all(flag in one for one in self._seq_flag)

    def pos2index(self: 'NAPosition', pos: int, first: str) -> int:
        if first == 'first':
            return self._seq_pos.index(pos)
        elif first == 'last':
            return len(self._seq_pos) - self._seq_pos[-1::-1].index(pos) - 1
        else:
            raise ValueError('invalid second parameter')

    def posrange2indexrange(
        self: 'NAPosition',
        pos_start: int,
        pos_end: int
    ) -> Tuple[int, int]:
        max_pos: int = self.max_pos()
        min_pos: int = self.min_pos()
        if max_pos < 0 or min_pos < 0:
            return 0, 0
        if pos_start > max_pos:
            idx_start = idx_end = self.pos2index(max_pos, 'last')
        elif pos_end < min_pos:
            idx_start = idx_end = self.pos2index(min_pos, 'first')
        else:
            pos_start = max(min_pos, pos_start)
            pos_end = min(max_pos, pos_end)
            for pos in range(pos_start, pos_end + 1):
                try:
                    idx_start = self.pos2index(pos, 'first')
                    break
                except ValueError:
                    continue
            for pos in range(pos_end, pos_start - 1, -1):
                try:
                    idx_end = self.pos2index(pos, 'last') + 1
                    break
                except ValueError:
                    continue
        return idx_start, idx_end

    def __getitem__(
        self: 'NAPosition',
        index: Union[int, slice]
    ) -> 'NAPosition':
        mytype = type(self)
        if isinstance(index, slice):
            return mytype(
                self._seq_text[index],
                self._seq_pos[index],
                self._seq_flag[index])
        else:
            return mytype(
                bytes(self._seq_text[index]),
                [self._seq_pos[index]],
                [self._seq_flag[index]])

    def __setitem__(
        self: 'NAPosition',
        index: Union[int, slice],
        value: Union['NAPosition', bytes, None]
    ) -> None:
        if isinstance(index, int):
            index = slice(index, index + 1 or None)
        if isinstance(value, NAPosition):
            seq_bytes = bytearray(self._seq_text)
            seq_bytes[index] = bytearray(value._seq_text)
            self._seq_text = bytes(seq_bytes)
            self._seq_pos[index] = value._seq_pos
            self._seq_flag[index] = value._seq_flag
        elif not value:
            seq_bytes = bytearray(self._seq_text)
            seq_bytes[index] = b''
            self._seq_text = bytes(seq_bytes)
            self._seq_pos[index] = []
            self._seq_flag[index] = []
        elif isinstance(value, bytes) and \
                all(na in GAP_CHARS for na in value):
            seq_bytes = bytearray(self._seq_text)
            seq_bytes[index] = bytearray(value)
            self._seq_text = bytes(seq_bytes)
            self._seq_pos[index] = [-1] * len(value)
            self._seq_flag[index] = [set() for _ in value]
        else:
            raise ValueError('invalid input')

    def __len__(self: 'NAPosition') -> int:
        return len(self._seq_text)

    def __contains__(self: 'NAPosition', value: bytes) -> bool:
        return value in self._seq_text

    def __str__(self: 'NAPosition') -> str:
        return str(self._seq_text, 'ASCII')

    def __bytes__(self: 'NAPosition') -> bytes:
        return self._seq_text

    def __repr__(self: 'NAPosition') -> str:
        return 'NAPosition({!r})'.format(self._seq_text)

    def __iter__(self: 'NAPosition') -> Generator['NAPosition', None, None]:
        na: int
        pos: int
        flag: Set[str]
        for na, pos, flag in zip(self._seq_text,
                                 self._seq_pos,
                                 self._seq_flag):
            yield NAPosition(bytes([na]), [pos], [flag])

    def __iadd__(self: 'NAPosition', other: 'NAPosition') -> None:
        self._seq_text += other._seq_text
        self._seq_pos += other._seq_pos
        self._seq_flag += other._seq_flag

    def __add__(self: 'NAPosition', other: 'NAPosition') -> 'NAPosition':
        seq_text = self._seq_text + other._seq_text
        seq_pos = self._seq_pos + other._seq_pos
        seq_flag = self._seq_flag + other._seq_flag
        return type(self)(seq_text, seq_pos, seq_flag)

    # @cython.ccall
    def count(
        self: 'NAPosition',
        sub: Union[int, bytes],
        start: Optional[int] = None,
        end: Optional[int] = None
    ) -> int:
        return self._seq_text.count(sub, start, end)

    # @cython.ccall
    def remove_gaps(self: 'NAPosition') -> 'NAPosition':
        return type(self).init_from_triplets([
            (na, pos, flag)
            for na, pos, flag in zip(self._seq_text,
                                     self._seq_pos,
                                     self._seq_flag)
            if na not in GAP_CHARS
        ])

    # @cython.ccall
    def is_gap(self: 'NAPosition') -> bool:
        if self._is_gap is None:
            self._is_gap = False
            if self._seq_text:
                gap: int
                total_gap: int = 0
                for gap in GAP_CHARS:
                    total_gap += self._seq_text.count(gap)
                self._is_gap = total_gap == len(self._seq_text)
        return self._is_gap

    @classmethod
    def join(
        cls: Type['NAPosition'],
        value: Iterable['NAPosition']
    ) -> 'NAPosition':
        seq_text: bytes = b''.join([one._seq_text for one in value])
        seq_pos: List[int] = [pos for one in value for pos in one._seq_pos]
        seq_flag: List[Set[str]] = [
            flag for one in value for flag in one._seq_flag
        ]
        return cls(seq_text, seq_pos, seq_flag)

    @classmethod
    def join_as_bytes(
        cls: Type['NAPosition'],
        value: Iterable['NAPosition']
    ) -> bytes:
        return b''.join([one._seq_text for one in value])

    @classmethod
    def join_as_str(
        cls: Type['NAPosition'],
        value: Iterable['NAPosition']
    ) -> str:
        seq_text: bytes = b''.join([one._seq_text for one in value])
        return str(seq_text, 'ASCII')

    @staticmethod
    def list_contains_any_gap(nas: Iterable['NAPosition']) -> bool:
        for na in nas:
            if na.is_gap():
                return True
        return False

    @staticmethod
    def list_contains_all_gap(nas: Iterable['NAPosition']) -> bool:
        for na in nas:
            if not na.is_gap():
                return False
        return True


NAPosOrList = TypeVar('NAPosOrList', NAPosition, List[NAPosition])
