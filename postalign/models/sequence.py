import re
from collections import namedtuple

from .modifier import ModifierLinkedList

GAP_CHARS = '.-'

SEQTYPES = ['NA']  # currently we don't support AA sequence

SKIP_VALIDATION = object()

ILLEGAL_PATTERNS = {
    'NA': re.compile(r'[^ACGTUWSMKRYBDHVN.-]')
}

BaseSequence = namedtuple(
    'BaseSequence', [
        'header',
        'description',
        'seqtext',
        'seqid',
        'seqtype',
        # absolute seqstart (>=0) according to original cmd input
        'abs_seqstart',
        'modifiers_'
    ]
)


class Sequence(BaseSequence):

    def __new__(
        cls, *,
        header,
        description,
        seqtext,
        seqid,
        seqtype,
        abs_seqstart,
        modifiers_=None,
        skip_invalid=True
    ):
        seqtext = seqtext.upper()
        illegal_pattern = ILLEGAL_PATTERNS[seqtype]
        invalids = illegal_pattern.findall(seqtext)
        if skip_invalid is not SKIP_VALIDATION:
            if invalids and skip_invalid:
                seqtext = illegal_pattern.sub('-', seqtext)
            elif invalids:
                raise ValueError(
                    'sequence {} contains invalid notation(s) ({})'
                    'while skip_invalid=False'
                    .format(header, ''.join(set(invalids)))
                )
        return super().__new__(
            cls,
            header,
            description,
            seqtext,
            seqid,
            seqtype,
            abs_seqstart,
            modifiers_)

    @property
    def modifiers(self):
        if self.modifiers_ is None:
            return ModifierLinkedList()
        else:
            return self.modifiers_

    def __getitem__(self, index):
        seqtext = self.seqtext[index]
        if isinstance(index, slice):
            replace_flag = True
            sliceval = index.indices(len(self.seqtext))
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

            return Sequence(
                header=self.header,
                description=self.description,
                seqtext=seqtext,
                seqid=self.seqid,
                seqtype=self.seqtype,
                modifiers_=modifiers,
                abs_seqstart=abs_seqstart,
                skip_invalid=SKIP_VALIDATION)
        return seqtext

    def __add__(self, other):
        if not isinstance(other, BaseSequence):
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
        return Sequence(
            header=self.header,
            seqtext=self.seqtext + other.seqtext,
            seqid=self.seqid,
            seqtype=self.seqtype,
            modifiers_=self.modifiers + other.modifiers,
            abs_seqstart=self.abs_seqstart,
            skip_invalid=SKIP_VALIDATION)

    def __iter__(self):
        yield from self.seqtext

    def __len__(self):
        return len(self.seqtext)

    def push_seqtext(self, seqtext, modtext, start_offset, **kw):
        """Modify seqtext and push modifier forward

        The difference between `push_seqtext` and `replace_seqtext` is similar
        to modern browsers' `history.pushstate` and `history.replacestate`,
        where the history is stored in attribute `modifiers`.
        """
        modifiers = self.modifiers.push(modtext, **kw)
        abs_seqstart = self.abs_seqstart + start_offset

        return Sequence(
            header=self.header,
            description=self.description,
            seqtext=seqtext,
            seqid=self.seqid,
            seqtype=self.seqtype,
            modifiers_=modifiers,
            abs_seqstart=abs_seqstart,
            skip_invalid=True)

    def replace_seqtext(self, seqtext, modtext, start_offset, **kw):
        """Modify seqtext and replace last modifier

        The difference between `push_seqtext` and `replace_seqtext` is similar
        to modern browsers' `history.pushstate` and `history.replacestate`,
        where the history is stored in attribute `modifiers`.
        """
        modifiers = self.modifiers.replace_last(modtext, **kw)
        abs_seqstart = self.abs_seqstart + start_offset

        return Sequence(
            header=self.header,
            description=self.description,
            seqtext=seqtext,
            seqid=self.seqid,
            seqtype=self.seqtype,
            modifiers_=modifiers,
            abs_seqstart=abs_seqstart,
            skip_invalid=True)

    @property
    def header_with_modifiers(self):
        modtext = str(self.modifiers)
        if modtext:
            return '{} MOD::{}'.format(self.header, self.modifiers)
        else:
            return self.header
