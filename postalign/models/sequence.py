import re
from typing import (
    Tuple,
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

GAP_CHARS = '.-'
GAP_PATTERN = re.compile(r'^[.-]+$')
ANY_GAP_PATTERN = re.compile(r'[.-]')

SEQTYPES = ['NA']  # currently we don't support AA sequence

SKIP_VALIDATION = object()

ILLEGAL_PATTERNS = {
    'NA': re.compile(r'[^ACGTUWSMKRYBDHVN.-]')
}

T = TypeVar('T', bound='Sequence')
Position = TypeVar('Position', NAPosition, AAPosition)


class Sequence(Generic[Position]):
    header: str
    description: str
    seqtext: Position
    seqid: int
    seqtype: str
    abs_seqstart: int
    modifiers_: Optional[ModifierLinkedList]

    def __init__(
        self: T, *,
        header: str,
        description: str,
        seqtext: Position,
        seqid: int,
        seqtype: str,
        abs_seqstart: int,
        modifiers_: Optional[ModifierLinkedList] = None,
        skip_invalid: Union[bool, object] = True
    ) -> None:
        if seqtype == 'AA' or isinstance(seqtext, AAPosition):
            raise NotImplementedError('Amino acid is not support yet')
        if seqtype != 'NA' and not isinstance(seqtext, NAPosition):
            raise ValueError(
                "seqtext must be an instance of "
                "NAPosition when seqtype is 'NA'"
            )
        if skip_invalid is not SKIP_VALIDATION:
            testseqtext: str = str(bytes(seqtext), 'U8')
            illegal_pattern: re.Pattern = ILLEGAL_PATTERNS[seqtype]
            invalids: List[str] = illegal_pattern.findall(testseqtext)
            if invalids and skip_invalid:
                seqtext = NAPosition.init_from_nastring(
                    bytes(illegal_pattern.sub('-', testseqtext), 'U8')
                )
            elif invalids:
                raise ValueError(
                    'sequence {} contains invalid notation(s) ({})'
                    'while skip_invalid=False'
                    .format(header, ''.join(set(invalids)))
                )
        self.header = header
        self.descrption = description
        self.seqtext = seqtext
        self.seqid = seqid
        self.seqtype = seqtype
        self.abs_seqstart = abs_seqstart
        self.modifiers_ = modifiers_

    @property
    def headerdesc(self: T) -> str:
        hd: str = self.header
        if self.description:
            hd += ' ' + self.description
        return hd

    @property
    def modifiers(self: T) -> ModifierLinkedList:
        if self.modifiers_ is None:
            return ModifierLinkedList()
        else:
            return self.modifiers_

    def __getitem__(self: T, index: Union[int, slice]) -> Union[Position, T]:
        seqtext: Position = self.seqtext[index]
        if isinstance(index, slice):
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
                    'step slicing is unsupported: {!r}'.format(index)
                )
            if replace_flag:
                modifiers = self.modifiers.replace_last(
                    modtext, slicetuples=slicetuples)
            else:
                modifiers = self.modifiers.push(
                    modtext, slicetuples=slicetuples)

            abs_seqstart = self.abs_seqstart + start
            for gap in GAP_CHARS:
                abs_seqstart -= self.seqtext[:start].count(gap)

            return type(self)(
                header=self.header,
                description=self.description,
                seqtext=seqtext,
                seqid=self.seqid,
                seqtype=self.seqtype,
                modifiers_=modifiers,
                abs_seqstart=abs_seqstart,
                skip_invalid=SKIP_VALIDATION)
        return seqtext

    def __add__(self: T, other: T) -> T:
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

    def __iter__(self: T) -> Generator[Position, None, None]:
        yield from self.seqtext

    def __len__(self: T) -> int:
        return len(self.seqtext)

    def push_seqtext(
        self: T,
        seqtext: Position,
        modtext: str,
        start_offset: int,
        **kw: Any
    ) -> T:
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
        self: T,
        seqtext: Position,
        modtext: str,
        start_offset: int,
        **kw: Any
    ) -> T:
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
    def header_with_modifiers(self: T) -> str:
        modtext: str = str(self.modifiers)
        if modtext:
            return '{} MOD::{}'.format(self.header, self.modifiers)
        else:
            return self.header


RefSeqPair = Tuple[Sequence, Sequence]
