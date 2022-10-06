from typing import Type, List, ByteString


class AAPosition:

    notation: int
    pos: int
    flag: int
    is_gap: bool

    @classmethod
    def init_gaps(cls: Type['AAPosition'], gaplen: int) -> List['AAPosition']:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @classmethod
    def init_from_bytes(
        cls: Type['AAPosition'],
        seq_text: ByteString
    ) -> List['AAPosition']:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def any_has_gap(aas: List['AAPosition']) -> bool:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def all_have_gap(aas: List['AAPosition']) -> bool:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def count_gaps(nas: List['AAPosition']) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def count_nongaps(nas: List['AAPosition']) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @classmethod
    def as_str(
        cls: Type['AAPosition'],
        nas: List['AAPosition']
    ) -> str:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def min_pos(nas: List['AAPosition']) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def max_pos(nas: List['AAPosition']) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def min_nongap_index(nas: List['AAPosition']) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def max_nongap_index(nas: List['AAPosition']) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')
