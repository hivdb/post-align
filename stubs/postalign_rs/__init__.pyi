"""Type stubs for the postalign_rs Rust extension (PyO3)."""

from __future__ import annotations

from typing import Generic, TypeVar, overload

_T = TypeVar('_T')

# ---------------------------------------------------------------------------
# NAPosition
# ---------------------------------------------------------------------------


class NAPosition(Generic[_T]):
    notation: int
    pos: int
    flag: int
    is_gap: bool
    payload: _T | None

    def __init__(
        self,
        notation: int,
        pos: int,
        flag: int,
        payload: _T | None = None,
    ) -> None: ...
    def __str__(self) -> str: ...
    def __bytes__(self) -> bytes: ...
    def __repr__(self) -> str: ...
    def __copy__(self) -> NAPosition[_T]: ...

    # -- Class / static methods operating on list[NAPosition] ---------------

    @classmethod
    def init_gaps(cls, gaplen: int) -> list[NAPosition[None]]: ...
    @classmethod
    def init_from_bytes(
        cls,
        seq_text: bytes,
        seq_payload: list[_T] | None = None,
    ) -> list[NAPosition[_T]]: ...
    @staticmethod
    def min_pos(nas: list[NAPosition[_T]]) -> int: ...
    @staticmethod
    def max_pos(nas: list[NAPosition[_T]]) -> int: ...
    @staticmethod
    def min_nongap_index(
        nas: list[NAPosition[_T]], start: int = -1, stop: int = -1
    ) -> int: ...
    @staticmethod
    def max_nongap_index(
        nas: list[NAPosition[_T]], start: int = -1, stop: int = -1
    ) -> int: ...
    @staticmethod
    def set_flag(nas: list[NAPosition[_T]], flag: int) -> None: ...
    @staticmethod
    def any_has_flag(nas: list[NAPosition[_T]], flag: int) -> bool: ...
    @staticmethod
    def all_have_flag(nas: list[NAPosition[_T]], flag: int) -> bool: ...
    @staticmethod
    def count_gaps(nas: list[NAPosition[_T]]) -> int: ...
    @staticmethod
    def count_nongaps(nas: list[NAPosition[_T]]) -> int: ...
    @staticmethod
    def any_has_gap(nas: list[NAPosition[_T]]) -> bool: ...
    @staticmethod
    def all_have_gap(nas: list[NAPosition[_T]]) -> bool: ...
    @staticmethod
    def remove_gaps(nas: list[NAPosition[_T]]) -> list[NAPosition[_T]]: ...
    @staticmethod
    def as_bytes(nas: list[NAPosition[_T]]) -> bytes: ...
    @classmethod
    def as_str(cls, nas: list[NAPosition[_T]]) -> str: ...
    @staticmethod
    def posrange2indexrange(
        nas: list[NAPosition[_T]],
        pos_start: int,
        pos_end: int,
        include_boundary_gaps: bool = False,
    ) -> tuple[int, int]: ...


# ---------------------------------------------------------------------------
# NAPositionListIter
# ---------------------------------------------------------------------------


class NAPositionListIter:
    def __iter__(self) -> NAPositionListIter: ...
    def __next__(self) -> NAPosition[None]: ...


# ---------------------------------------------------------------------------
# NAPositionList
# ---------------------------------------------------------------------------


class NAPositionList:
    def __init__(self) -> None: ...
    def __len__(self) -> int: ...
    def __bool__(self) -> bool: ...
    def __repr__(self) -> str: ...
    @overload
    def __getitem__(self, index: int) -> NAPosition[None]: ...
    @overload
    def __getitem__(self, index: slice) -> NAPositionList: ...
    def __iter__(self) -> NAPositionListIter: ...
    def __add__(self, other: NAPositionList) -> NAPositionList: ...

    # -- Constructors -------------------------------------------------------

    @classmethod
    def from_bytes(cls, seq_text: bytes) -> NAPositionList: ...
    @classmethod
    def from_arrays(
        cls,
        notations: list[int],
        positions: list[int],
        flags: list[int],
    ) -> NAPositionList: ...
    @classmethod
    def init_gaps(cls, gaplen: int) -> NAPositionList: ...
    @classmethod
    def from_list(cls, nas: list[NAPosition[None]]) -> NAPositionList: ...
    def to_list(self) -> list[NAPosition[None]]: ...

    # -- Bulk read ----------------------------------------------------------

    def count_gaps(self) -> int: ...
    def count_nongaps(self) -> int: ...
    def any_has_gap(self) -> bool: ...
    def all_have_gap(self) -> bool: ...
    def any_has_flag(self, flag: int) -> bool: ...
    def all_have_flag(self, flag: int) -> bool: ...
    def as_bytes(self) -> bytes: ...
    def as_str(self) -> str: ...

    # -- Bulk write ---------------------------------------------------------

    def set_flag(self, flag: int) -> None: ...

    # -- Index operations ---------------------------------------------------

    def min_pos(self) -> int: ...
    def max_pos(self) -> int: ...
    def min_nongap_index(self, start: int = -1, stop: int = -1) -> int: ...
    def max_nongap_index(self, start: int = -1, stop: int = -1) -> int: ...
    def remove_gaps(self) -> NAPositionList: ...
    def posrange2indexrange(
        self,
        pos_start: int,
        pos_end: int,
        include_boundary_gaps: bool = False,
    ) -> tuple[int, int]: ...


# ---------------------------------------------------------------------------
# AAPosition
# ---------------------------------------------------------------------------


class AAPosition(Generic[_T]):
    notation: int
    pos: int
    flag: int
    is_gap: bool
    payload: _T | None

    def __init__(
        self,
        notation: int,
        pos: int,
        flag: int,
        payload: _T | None = None,
    ) -> None: ...
    def __copy__(self) -> AAPosition[_T]: ...

    @classmethod
    def init_gaps(cls, gaplen: int) -> list[AAPosition[None]]: ...
    @classmethod
    def init_from_bytes(
        cls,
        seq_text: bytes,
        seq_payload: list[_T] | None = None,
    ) -> list[AAPosition[_T]]: ...


# ---------------------------------------------------------------------------
# Module-level functions
# ---------------------------------------------------------------------------


def enumerate_seq_pos(seq_text: bytes) -> list[int]: ...
def realign_gaps(
    ref_notations: list[int],
    ref_positions: list[int],
    ref_flags: list[int],
    seq_notations: list[int],
    seq_positions: list[int],
    seq_flags: list[int],
    min_gap_distance: int,
    window_size: int,
    gap_placement_score: dict[int, dict[tuple[int, int], int]],
    is_seq_start: bool,
    is_seq_end: bool,
) -> tuple[list[int], list[int], list[int], list[int], list[int], list[int]]: ...
def realign_gaps_optimized(
    ref_notations: list[int],
    ref_positions: list[int],
    ref_flags: list[int],
    seq_notations: list[int],
    seq_positions: list[int],
    seq_flags: list[int],
    min_gap_distance: int,
    window_size: int,
    gap_placement_score: dict[int, dict[tuple[int, int], int]],
    is_seq_start: bool,
    is_seq_end: bool,
) -> tuple[list[int], list[int], list[int], list[int], list[int], list[int]]: ...
def codon_align_full(
    ref_notations: list[int],
    ref_positions: list[int],
    ref_flags: list[int],
    seq_notations: list[int],
    seq_positions: list[int],
    seq_flags: list[int],
    min_gap_distance: int,
    window_size: int,
    gap_placement_score: dict[int, dict[tuple[int, int], int]],
    ref_start: int,
    ref_end: int,
) -> tuple[int, int, list[int], list[int], list[int], list[int], list[int], list[int]] | None: ...
def codon_align_full_v2(
    ref_nas: NAPositionList,
    seq_nas: NAPositionList,
    min_gap_distance: int,
    window_size: int,
    gap_placement_score: dict[int, dict[tuple[int, int], int]],
    ref_start: int,
    ref_end: int,
) -> tuple[int, int, NAPositionList, NAPositionList] | None: ...
def codon_align_batch(
    items: list[
        tuple[list[int], list[int], list[int], list[int], list[int], list[int], int, int]
    ],
    min_gap_distance: int,
    window_size: int,
    gap_placement_score: dict[int, dict[tuple[int, int], int]],
) -> list[
    tuple[int, int, list[int], list[int], list[int], list[int], list[int], list[int]] | None
]: ...
