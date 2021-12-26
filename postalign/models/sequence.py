import re
from typing import (
    Tuple,
    Type,
    TypeVar,
    Generic,
    Optional,
    List,
    Union,
    Any,
    Generator
)

from .na_position import NAPosition
from .aa_position import AAPosition

from .modifier import ModifierLinkedList
from ._sequence import sanitize_sequence, SKIP_VALIDATION

GAP_CHARS = '.-'
GAP_PATTERN = re.compile(r'^[.-]+$')
ANY_GAP_PATTERN = re.compile(r'[.-]')

Position = TypeVar('Position', NAPosition, AAPosition)


class Sequence(Generic[Position]):
    header: str
    description: str
    seqtext: List[Position]
    seqid: int
    seqtype: Type[Position]
    abs_seqstart: int
    modifiers_: Optional[ModifierLinkedList]

    def __init__(
        self: 'Sequence', *,
        header: str,
        description: str,
        seqtext: List[Position],
        seqid: int,
        seqtype: Type[Position],
        abs_seqstart: int,
        modifiers_: Optional[ModifierLinkedList] = None,
        skip_invalid: Union[bool, object] = True
    ) -> None:
        one: Position

        if skip_invalid is not SKIP_VALIDATION:
            seqtext = sanitize_sequence(seqtext, seqtype, header, skip_invalid)

        self.header = header
        self.description = description
        self.seqtext = seqtext
        self.seqid = seqid
        self.seqtype = seqtype
        self.abs_seqstart = abs_seqstart
        self.modifiers_ = modifiers_

    @property
    def seqtext_as_str(self: 'Sequence') -> str:
        seqtype: Type[Position] = self.seqtype
        return seqtype.as_str(self.seqtext)

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

    def __getitem__(
        self: 'Sequence',
        index: Union[int, slice]
    ) -> Union[Position, 'Sequence']:
        seqtext: List[Position] = self.seqtext
        if isinstance(index, slice):
            seqtext = seqtext[index]
            seqtype: Type[Position] = self.seqtype
            start: int
            end: int
            slice_remain_len: int
            accum_len: int
            slicetuples: List[Tuple[int, int]]
            prevslicetuples: List[Tuple[int, int]]
            modtext: str
            modifiers: ModifierLinkedList
            replace_flag: bool = True
            sliceval: Tuple[int, int, int] = index.indices(len(self.seqtext))
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
                    modtext = 'join({})'.format(
                        ','.join('{}..{}'.format(*stuple)
                                 for stuple in slicetuples)
                    )
            else:
                raise ValueError(
                    'step slicing is not supported: {!r}'.format(index)
                )
            if replace_flag:
                modifiers = self.modifiers.replace_last(
                    modtext, slicetuples=slicetuples)
            else:
                modifiers = self.modifiers.push(
                    modtext, slicetuples=slicetuples)

            abs_seqstart = self.abs_seqstart + start
            for gap in GAP_CHARS:
                abs_seqstart -= seqtype.count_gaps(self.seqtext[:start])

            return Sequence(
                header=self.header,
                description=self.description,
                seqtext=seqtext,
                seqid=self.seqid,
                seqtype=self.seqtype,
                modifiers_=modifiers,
                abs_seqstart=abs_seqstart,
                skip_invalid=SKIP_VALIDATION)
        else:
            return seqtext[index]

    def __add__(self: 'Sequence', other: 'Sequence') -> 'Sequence':
        if not isinstance(other, Sequence):
            raise TypeError(
                'unsupported operand type(s) for +: {!r} and {!r}'
                .format(self.__class__.__name__, other.__class__.__name__)
            )
        if self.seqid != other.seqid:
            raise ValueError(
                'concat two sequences with different seqid is disallowed: '
                '{!r} and {!r}'
                .format(self.header, other.header)
            )
        return type(self)(
            header=self.header,
            description=self.description,
            seqtext=self.seqtext + other.seqtext,
            seqid=self.seqid,
            seqtype=self.seqtype,
            modifiers_=self.modifiers + other.modifiers,
            abs_seqstart=self.abs_seqstart,
            skip_invalid=SKIP_VALIDATION)

    def __iter__(self: 'Sequence') -> Generator[Position, None, None]:
        yield from self.seqtext

    def __len__(self: 'Sequence') -> int:
        return len(self.seqtext)

    def push_seqtext(
        self: 'Sequence',
        seqtext: List[Position],
        modtext: str,
        start_offset: int,
        **kw: Any
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
            skip_invalid=True)

    def replace_seqtext(
        self: 'Sequence',
        seqtext: List[Position],
        modtext: str,
        start_offset: int,
        **kw: Any
    ) -> 'Sequence':
        """Modify seqtext and replace last modifier

        The difference between `push_seqtext` and `replace_seqtext` is similar
        to modern browsers' `history.pushstate` and `history.replacestate`,
        where the history is stored in attribute `modifiers`.
        """
        modifiers: ModifierLinkedList = \
            self.modifiers.replace_last(modtext, **kw)
        abs_seqstart: int = self.abs_seqstart + start_offset

        return type(self)(
            header=self.header,
            description=self.description,
            seqtext=seqtext,
            seqid=self.seqid,
            seqtype=self.seqtype,
            modifiers_=modifiers,
            abs_seqstart=abs_seqstart,
            skip_invalid=True)

    @property
    def header_with_modifiers(self: 'Sequence') -> str:
        modtext: str = str(self.modifiers)
        if modtext:
            return '{} MOD::{}'.format(self.header, self.modifiers)
        else:
            return self.header


RefSeqPair = Tuple[Sequence, Sequence]
