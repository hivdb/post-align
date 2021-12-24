import cython  # type: ignore
from typing import (
    List,
    Set,
    Optional,
    Tuple,
    Union,
    Iterable,
    Generator,
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


T = TypeVar('T', bound='NAPosition')


class NAPosition:

    @classmethod
    def init_empty(cls: Type[T]) -> T:
        return cls(b'', [], [])

    @classmethod
    def init_gaps(cls: Type[T], gaplen: int) -> T:
        seq_text: bytes = b'-' * gaplen
        seq_pos: List[int] = [-1] * gaplen
        seq_flag: List[Set[str]] = [set() for _ in seq_pos]
        return cls(seq_text, seq_pos, seq_flag)

    @classmethod
    def init_from_triplets(
        cls: Type[T],
        triplets: List[Tuple[int, int, Set[str]]]
    ) -> T:
        seq_text: bytearray = bytearray()
        seq_pos: List[int] = []
        seq_flag: List[Set[str]] = []
        for na, pos, flag in triplets:
            seq_text.append(na)
            seq_pos.append(pos)
            seq_flag.append(flag)
        return cls(seq_text, seq_pos, seq_flag)

    @classmethod
    def init_from_nastring(cls: Type[T], seq_text: bytes) -> T:
        seq_text = seq_text.upper()
        seq_pos: List[int] = enumerate_seq_pos(seq_text)
        seq_flag: List[Set[str]] = [set() for _ in seq_pos]
        return cls(seq_text, seq_pos, seq_flag)

    _seq_text: bytes
    _seq_pos: List[int]
    _seq_flag: List[Set[str]]
    _min_pos: Optional[int]
    _max_pos: Optional[int]
    _is_gap: Optional[bool]
    _is_single_gap: Optional[bool]

    def __init__(
        self: T,
        seq_text: bytes,
        seq_pos: List[int],
        seq_flag: List[Set[str]]
    ) -> None:
        self._seq_text = seq_text
        self._seq_pos = seq_pos
        self._seq_flag = seq_flag
        self._min_pos = None
        self._max_pos = None
        if len(seq_text) == 1:
            self._is_gap = seq_text in GAP_CHARS
            self.is_single_gap = self._is_gap
        else:
            self._is_gap = None
            self.is_single_gap = False

    @property
    def min_pos(self: T) -> Optional[int]:
        if self._min_pos is None:
            self._min_pos = next(
                (pos for pos in self._seq_pos if pos > 0), None
            )
        return self._min_pos

    @property
    def max_pos(self: T) -> Optional[int]:
        if self._max_pos is None:
            self._max_pos = next(
                (pos for pos in reversed(self._seq_pos) if pos > 0), None
            )
        return self._max_pos

    @property
    def empty(self: T) -> bool:
        return self.min_pos is None

    @cython.ccall
    def set_flag(self: T, flag: str) -> None:
        one: Set[str]
        for one in self._seq_flag:
            one.add(flag)

    @cython.ccall
    def any_has_flag(self: T, flag: str) -> bool:
        one: Set[str]
        for one in self._seq_flag:
            if flag in one:
                return True
        return False

    @cython.ccall
    def all_have_flag(self: T, flag: str) -> bool:
        one: Set[str]
        return all(flag in one for one in self._seq_flag)

    @cython.ccall
    def pos2index(self: T, pos: int, first: str) -> int:
        if first == 'first':
            return self._seq_pos.index(pos)
        elif first == 'last':
            return len(self._seq_pos) - self._seq_pos[-1::-1].index(pos) - 1
        else:
            raise ValueError('invalid second parameter')

    @cython.ccall
    def posrange2indexrange(
        self: T,
        pos_start: int,
        pos_end: int
    ) -> Tuple[int, int]:
        max_pos: Optional[int] = self.max_pos
        min_pos: Optional[int] = self.min_pos
        if max_pos is None or min_pos is None:
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

    def __getitem__(self: T, index: Union[int, slice]) -> T:
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
        self: T,
        index: Union[int, slice],
        value: Union[T, bytes, None]
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

    def __len__(self: T) -> int:
        return len(self._seq_text)

    def __contains__(self: T, value: bytes) -> bool:
        return value in self._seq_text

    def __str__(self: T) -> str:
        raise NotImplementedError(
            'str() for NAPosition is deliberately not supported, '
            'use bytes() instead'
        )

    def __bytes__(self: T) -> bytes:
        return self._seq_text

    def __repr__(self: T) -> str:
        return 'NAPosition({!r})'.format(self._seq_text)

    def __iter__(self: T) -> Generator[T, None, None]:
        na: int
        pos: int
        flag: Set[str]
        mytype = type(self)
        for na, pos, flag in zip(self._seq_text,
                                 self._seq_pos,
                                 self._seq_flag):
            yield mytype(bytes(na), [pos], [flag])

    def __iadd__(self: T, other: T) -> None:
        self._seq_text += other._seq_text
        self._seq_pos += other._seq_pos
        self._seq_flag += other._seq_flag

    def __add__(self: T, other: T) -> T:
        seq_text = self._seq_text + other._seq_text
        seq_pos = self._seq_pos + other._seq_pos
        seq_flag = self._seq_flag + other._seq_flag
        return type(self)(seq_text, seq_pos, seq_flag)

    @cython.ccall
    def count(
        self: T,
        sub: Union[int, bytes],
        start: Optional[int] = None,
        end: Optional[int] = None
    ) -> int:
        return self._seq_text.count(sub, start, end)

    @cython.ccall
    def remove_gaps(self: T) -> T:
        return type(self).init_from_triplets([
            (na, pos, flag)
            for na, pos, flag in zip(self._seq_text,
                                     self._seq_pos,
                                     self._seq_flag)
            if na not in GAP_CHARS
        ])

    @cython.ccall
    def is_gap(self: T) -> bool:
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
    def join(cls: Type[T], value: Iterable[T]) -> T:
        seq_text: bytes = b''.join([one._seq_text for one in value])
        seq_pos: List[int] = [pos for one in value for pos in one._seq_pos]
        seq_flag: List[Set[str]] = [
            flag for one in value for flag in one._seq_flag
        ]
        return cls(seq_text, seq_pos, seq_flag)

    @staticmethod
    def list_contains_any_gap(nas: Iterable[T]) -> bool:
        for na in nas:
            if na.is_gap():
                return True
        return False

    @staticmethod
    def list_contains_all_gap(nas: Iterable[T]) -> bool:
        for na in nas:
            if not na.is_gap():
                return False
        return True


NAPosOrList = TypeVar('NAPosOrList', NAPosition, List[NAPosition])