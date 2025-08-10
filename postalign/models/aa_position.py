import cython  # type: ignore
from typing import Any
from collections.abc import ByteString
from .position_flag import PositionFlag


@cython.cclass
class AAPosition:

    notation: int = cython.declare(cython.int, visibility='public')
    pos: int = cython.declare(cython.int, visibility='public')
    flag: PositionFlag = cython.declare(cython.int, visibility='public')
    is_gap: bool = cython.declare(cython.bint, visibility='public')

    payload: Any = cython.declare(object, visibility='public')

    def __copy__(self: 'AAPosition') -> 'AAPosition':
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @classmethod
    def init_gaps(cls: type['AAPosition'], gaplen: int) -> list['AAPosition']:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @classmethod
    def init_from_bytes(
        cls: type['AAPosition'],
        seq_text: ByteString,
        seq_payload: list | None = None
    ) -> list['AAPosition']:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def any_has_gap(aas: list['AAPosition']) -> bool:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def all_have_gap(aas: list['AAPosition']) -> bool:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def count_gaps(nas: list['AAPosition']) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def count_nongaps(nas: list['AAPosition']) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @classmethod
    def as_str(
        cls: type['AAPosition'],
        nas: list['AAPosition']
    ) -> str:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def min_pos(nas: list['AAPosition']) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def max_pos(nas: list['AAPosition']) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def min_nongap_index(nas: list['AAPosition']) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def max_nongap_index(nas: list['AAPosition']) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def set_flag(nas: list['AAPosition'], flag: PositionFlag) -> None:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def any_has_flag(nas: list['AAPosition'], flag: PositionFlag) -> bool:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def all_have_flag(nas: list['AAPosition'], flag: PositionFlag) -> bool:
        raise NotImplementedError('Amino acid sequence is not yet supported')
