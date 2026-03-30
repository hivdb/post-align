"""Type stubs for the postalign_rs Rust extension (PyO3)."""

from __future__ import annotations

from typing import overload

# ---------------------------------------------------------------------------
# NAPosition
# ---------------------------------------------------------------------------


class NAPosition:
    notation: int
    pos: int
    flag: int
    is_gap: bool
    payload: object | None

    def __init__(
        self,
        notation: int,
        pos: int,
        flag: int,
        payload: object | None = None,
    ) -> None: ...
    def __str__(self) -> str: ...
    def __bytes__(self) -> bytes: ...
    def __repr__(self) -> str: ...
    def __copy__(self) -> NAPosition: ...

    # -- Class / static methods operating on list[NAPosition] ---------------

    @classmethod
    def init_gaps(cls, gaplen: int) -> list[NAPosition]: ...
    @classmethod
    def init_from_bytes(
        cls,
        seq_text: bytes,
        seq_payload: list[object] | None = None,
    ) -> list[NAPosition]: ...
    @staticmethod
    def min_pos(nas: list[NAPosition]) -> int: ...
    @staticmethod
    def max_pos(nas: list[NAPosition]) -> int: ...
    @staticmethod
    def min_nongap_index(
        nas: list[NAPosition], start: int = -1, stop: int = -1
    ) -> int: ...
    @staticmethod
    def max_nongap_index(
        nas: list[NAPosition], start: int = -1, stop: int = -1
    ) -> int: ...
    @staticmethod
    def set_flag(nas: list[NAPosition], flag: int) -> None: ...
    @staticmethod
    def any_has_flag(nas: list[NAPosition], flag: int) -> bool: ...
    @staticmethod
    def all_have_flag(nas: list[NAPosition], flag: int) -> bool: ...
    @staticmethod
    def count_gaps(nas: list[NAPosition]) -> int: ...
    @staticmethod
    def count_nongaps(nas: list[NAPosition]) -> int: ...
    @staticmethod
    def any_has_gap(nas: list[NAPosition]) -> bool: ...
    @staticmethod
    def all_have_gap(nas: list[NAPosition]) -> bool: ...
    @staticmethod
    def remove_gaps(nas: list[NAPosition]) -> list[NAPosition]: ...
    @staticmethod
    def as_bytes(nas: list[NAPosition]) -> bytes: ...
    @classmethod
    def as_str(cls, nas: list[NAPosition]) -> str: ...
    @staticmethod
    def posrange2indexrange(
        nas: list[NAPosition],
        pos_start: int,
        pos_end: int,
        include_boundary_gaps: bool = False,
    ) -> tuple[int, int]: ...


# ---------------------------------------------------------------------------
# NAPositionListIter
# ---------------------------------------------------------------------------


class NAPositionListIter:
    def __iter__(self) -> NAPositionListIter: ...
    def __next__(self) -> NAPosition: ...


# ---------------------------------------------------------------------------
# NAPositionList
# ---------------------------------------------------------------------------


class NAPositionList:
    def __init__(self) -> None: ...
    def __len__(self) -> int: ...
    def __bool__(self) -> bool: ...
    def __repr__(self) -> str: ...
    @overload
    def __getitem__(self, index: int) -> NAPosition: ...
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
    def from_list(cls, nas: list[NAPosition]) -> NAPositionList: ...
    def to_list(self) -> list[NAPosition]: ...

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


class AAPosition:
    notation: int
    pos: int
    flag: int
    is_gap: bool
    payload: object | None

    def __init__(
        self,
        notation: int,
        pos: int,
        flag: int,
        payload: object | None = None,
    ) -> None: ...
    def __copy__(self) -> AAPosition: ...

    @classmethod
    def init_gaps(cls, gaplen: int) -> list[AAPosition]: ...
    @classmethod
    def init_from_bytes(
        cls,
        seq_text: bytes,
        seq_payload: list[object] | None = None,
    ) -> list[AAPosition]: ...


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
