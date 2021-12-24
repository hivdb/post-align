from typing import Type, TypeVar, Union, Iterable, Generator

T = TypeVar('T', bound='AAPosition')


class AAPosition:

    @classmethod
    def init_gaps(cls: Type[T], gaplen: int) -> T:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def list_contains_any_gap(aas: Iterable[T]) -> bool:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    @staticmethod
    def list_contains_all_gap(aas: Iterable[T]) -> bool:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    def __getitem__(self: T, index: Union[int, slice]) -> T:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    def __add__(self: T, other: T) -> T:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    def __len__(self: T) -> int:
        raise NotImplementedError('Amino acid sequence is not yet supported')

    def __iter__(self: T) -> Generator[T, None, None]:
        raise NotImplementedError('Amino acid sequence is not yet supported')
