from typing import Callable, Iterable, Generic, TypeVar
from ..models.sequence import RefSeqPair


ReturnType = TypeVar('ReturnType', Iterable[RefSeqPair], Iterable[str])


class Processor(Generic[ReturnType]):
    command_name: str
    is_output_command: bool
    _processor: Callable[
        [Iterable[RefSeqPair]],
        ReturnType
    ]

    def __init__(
        self,
        command_name: str,
        is_output_command: bool,
        processor: Callable[
            [Iterable[RefSeqPair]],
            ReturnType
        ]
    ) -> None:
        self.command_name = command_name
        self.is_output_command = is_output_command
        # XXX: use setattr to avoid mypy warning:
        # https://github.com/python/mypy/issues/2427
        setattr(self, '_processor', processor)

    def __call__(self, iterator: Iterable[RefSeqPair]) -> ReturnType:
        return self._processor(iterator)  # type: ignore


def output_processor(command_name: str) -> Callable[
    [Callable[[Iterable[RefSeqPair]], Iterable[str]]],
    Processor[Iterable[str]]
]:

    def wrapper(processor: Callable[
        [Iterable[RefSeqPair]], Iterable[str]
    ]) -> Processor[Iterable[str]]:
        return Processor(command_name, True, processor)

    return wrapper


def intermediate_processor(command_name: str) -> Callable[
    [Callable[[Iterable[RefSeqPair]], Iterable[RefSeqPair]]],
    Processor[Iterable[RefSeqPair]]
]:

    def wrapper(processor: Callable[
        [Iterable[RefSeqPair]], Iterable[RefSeqPair]
    ]) -> Processor[Iterable[RefSeqPair]]:
        return Processor(command_name, False, processor)

    return wrapper
