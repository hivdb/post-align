import re
from collections.abc import Generator
from typing import Generic, TypeVar, Union

from ._sequence import SKIP_VALIDATION, _SkipValidationSentinel, sanitize_sequence
from .aa_position import AAPosition
from .modifier import ModifierLinkedList
from .na_position import NAPosition, NAPositionList

GAP_CHARS = '.-'
GAP_PATTERN = re.compile(r'^[.-]+$')
ANY_GAP_PATTERN = re.compile(r'[.-]')

Position = TypeVar('Position', NAPosition, AAPosition)


class Sequence(Generic[Position]):  # noqa: UP046
    header: str
    description: str
    seqtext: list[Position] | NAPositionList
    seqid: int
    seqtype: type[Position]
    abs_seqstart: int
    modifiers_: ModifierLinkedList | None

    def __init__(
        self: 'Sequence',
        *,
        header: str,
        description: str,
        seqtext: list[Position] | NAPositionList,
        seqid: int,
        seqtype: type[Position],
        abs_seqstart: int,
        modifiers_: ModifierLinkedList | None = None,
        skip_invalid: bool | _SkipValidationSentinel = True,
    ) -> None:
        if skip_invalid is not SKIP_VALIDATION:
            seqtext = sanitize_sequence(seqtext, seqtype, header, skip_invalid)

        if isinstance(seqtext, list) and seqtype is NAPosition:
            seqtext = NAPositionList.from_list(seqtext)  # type: ignore[arg-type]

        self.header = header
        self.description = description
        self.seqtext = seqtext
        self.seqid = seqid
        self.seqtype = seqtype
        self.abs_seqstart = abs_seqstart
        self.modifiers_ = modifiers_

    @property
    def seqtext_as_str(self: 'Sequence') -> str:
        st = self.seqtext
        if isinstance(st, NAPositionList):
            return st.as_str()
        raise NotImplementedError('AA seqtext_as_str is not yet supported')

    @property
    def headerdesc(self: 'Sequence') -> str:
        hd: str = self.header
        if self.description:
            hd += ' ' + self.description
        return hd

    @property
    def modifiers(self: 'Sequence') -> ModifierLinkedList:
        if self.modifiers_ is None:
            return ModifierLinkedList()
        else:
            return self.modifiers_

    def __getitem__(self: 'Sequence', index: int | slice) -> Union[Position, 'Sequence']:
        seqtext = self.seqtext
        if isinstance(index, slice):
            seqtext = seqtext[index]
            start: int
            end: int
            slice_remain_len: int
            accum_len: int
            slicetuples: list[tuple[int, int]]
            prevslicetuples: list[tuple[int, int]]
            modtext: str
            modifiers: ModifierLinkedList
            replace_flag: bool = True
            sliceval: tuple[int, int, int] = index.indices(len(self.seqtext))
            if sliceval[2] == 1:
                prevslicetuples = self.modifiers.last_modifier.slicetuples
                if not prevslicetuples:
                    replace_flag = False
                    prevslicetuples = [(0, len(self.seqtext))]
                start, end = sliceval[:2]
                slice_remain_len = end - start
                accum_len = 0
                slicetuples = []
                for frag_start, frag_end in prevslicetuples:
                    frag_len = frag_end - frag_start
                    accum_len += frag_len
                    if start < accum_len:
                        frag_start += start
                    if slice_remain_len <= frag_len:
                        frag_end = frag_start + slice_remain_len
                    slice_remain_len -= frag_len
                    slicetuples.append((frag_start, frag_end))
                if len(slicetuples) == 1:
                    modtext = 'slice({},{})'.format(*slicetuples[0])
                else:
                    modtext = 'join({})'.format(','.join('{}..{}'.format(*stuple) for stuple in slicetuples))
            else:
                raise ValueError(f'step slicing is not supported: {index!r}')
            if replace_flag:
                modifiers = self.modifiers.replace_last(modtext, slicetuples=slicetuples)
            else:
                modifiers = self.modifiers.push(modtext, slicetuples=slicetuples)

            abs_seqstart = self.abs_seqstart + start
            st_prefix = self.seqtext[:start]
            if isinstance(st_prefix, NAPositionList):
                abs_seqstart -= st_prefix.count_gaps()
            # else: AA path not yet supported

            return Sequence(
                header=self.header,
                description=self.description,
                seqtext=seqtext,
                seqid=self.seqid,
                seqtype=self.seqtype,
                modifiers_=modifiers,
                abs_seqstart=abs_seqstart,
                skip_invalid=SKIP_VALIDATION,
            )
        else:
            # seqtext is list[Position] | NAPositionList; both return NAPosition|AAPosition for int index.
            # mypy can't narrow the TypeVar union per-branch, so we suppress.
            return seqtext[index]  # type: ignore[return-value]

    def __add__(self: 'Sequence', other: 'Sequence') -> 'Sequence':
        if not isinstance(other, Sequence):
            raise TypeError(
                f'unsupported operand type(s) for +: {self.__class__.__name__!r} and {other.__class__.__name__!r}'
            )
        if self.seqid != other.seqid:
            raise ValueError(
                f'concat two sequences with different seqid is disallowed: {self.header!r} and {other.header!r}'
            )
        self_st = self.seqtext
        other_st = other.seqtext
        if isinstance(self_st, NAPositionList) and isinstance(other_st, NAPositionList):
            combined: list[Position] | NAPositionList = self_st + other_st
        elif isinstance(self_st, list) and isinstance(other_st, list):
            combined = self_st + other_st
        else:
            raise TypeError('cannot concatenate different seqtext types')
        return type(self)(
            header=self.header,
            description=self.description,
            seqtext=combined,
            seqid=self.seqid,
            seqtype=self.seqtype,
            modifiers_=self.modifiers + other.modifiers,
            abs_seqstart=self.abs_seqstart,
            skip_invalid=SKIP_VALIDATION,
        )

    def __iter__(self: 'Sequence') -> Generator[Position]:
        yield from self.seqtext

    def __len__(self: 'Sequence') -> int:
        return len(self.seqtext)

    def push_seqtext(
        self: 'Sequence',
        seqtext: list[Position] | NAPositionList,
        modtext: str,
        start_offset: int,
        **kw: list[tuple[int, int]],
    ) -> 'Sequence':
        """Modify seqtext and push modifier forward

        The difference between `push_seqtext` and `replace_seqtext` is similar
        to modern browsers' `history.pushstate` and `history.replacestate`,
        where the history is stored in attribute `modifiers`.
        """
        modifiers: ModifierLinkedList = self.modifiers.push(modtext, **kw)
        abs_seqstart: int = self.abs_seqstart + start_offset

        return type(self)(
            header=self.header,
            description=self.description,
            seqtext=seqtext,
            seqid=self.seqid,
            seqtype=self.seqtype,
            modifiers_=modifiers,
            abs_seqstart=abs_seqstart,
            skip_invalid=True,
        )

    def replace_seqtext(
        self: 'Sequence',
        seqtext: list[Position] | NAPositionList,
        modtext: str,
        start_offset: int,
        **kw: list[tuple[int, int]],
    ) -> 'Sequence':
        """Modify seqtext and replace last modifier

        The difference between `push_seqtext` and `replace_seqtext` is similar
        to modern browsers' `history.pushstate` and `history.replacestate`,
        where the history is stored in attribute `modifiers`.
        """
        modifiers: ModifierLinkedList = self.modifiers.replace_last(modtext, **kw)
        abs_seqstart: int = self.abs_seqstart + start_offset

        return type(self)(
            header=self.header,
            description=self.description,
            seqtext=seqtext,
            seqid=self.seqid,
            seqtype=self.seqtype,
            modifiers_=modifiers,
            abs_seqstart=abs_seqstart,
            skip_invalid=True,
        )

    @property
    def header_with_modifiers(self: 'Sequence') -> str:
        modtext: str = str(self.modifiers)
        if modtext:
            return f'{self.header} MOD::{self.modifiers}'
        else:
            return self.header


RefSeqPair = tuple[Sequence, Sequence]
