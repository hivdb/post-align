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


def enumerate_seq_pos(seq_text):
    offset = 1
    seq_pos = []
    for na in seq_text:
        if na in GAP_CHARS:
            seq_pos.append(-1)
        else:
            seq_pos.append(offset)
            offset += 1
    return seq_pos


class PositionalSeqStr:

    def __init__(self, seq_text, seq_pos=None, seq_flag=None):
        if not seq_text:
            self._seq_text = ''
            self._seq_pos = []
            self._seq_flag = []
        elif isinstance(seq_text, PositionalSeqStr):
            self._seq_text = seq_text._seq_text
            self._seq_pos = seq_text._seq_pos
            self._seq_flag = seq_text._seq_flag
        elif isinstance(seq_text, list):
            if isinstance(seq_text[0], tuple) and len(seq_text[0]) == 3:
                (self._seq_text,
                 self._seq_pos,
                 self._seq_flag) = list(zip(*seq_text))
            elif isinstance(seq_text[0], PositionalSeqStr):
                temp = PositionalSeqStr.join(seq_text)
                self._seq_text = temp._seq_text
                self._seq_pos = temp._seq_pos
                self._seq_flag = temp._seq_flag
            else:
                raise ValueError('invalid input')
        elif seq_pos is None and seq_flag is None:
            # used by utils.cigar
            self._seq_text = seq_text.upper()
            self._seq_pos = enumerate_seq_pos(seq_text)
            self._seq_flag = [set() for _ in self._seq_pos]
        elif isinstance(seq_pos, int) and isinstance(seq_flag, set):
            # used by __getitem__
            self._seq_text = seq_text.upper()
            self._seq_pos = [seq_pos]
            self._seq_flag = [seq_flag]
        elif len(seq_text) == len(seq_pos) == len(seq_flag):
            self._seq_text = seq_text.upper()
            self._seq_pos = seq_pos
            self._seq_flag = seq_flag
        else:
            raise ValueError('invalid input')
        non_gap_pos = [pos for pos in self._seq_pos if pos > 0]
        self.min_pos = self.max_pos = None
        if non_gap_pos:
            self.empty = False
            self.min_pos = min(non_gap_pos)
            self.max_pos = max(non_gap_pos)
        else:
            self.empty = True

    def set_flag(self, flag):
        for one in self._seq_flag:
            one.add(flag)

    def any_has_flag(self, flag):
        return any(flag in one for one in self._seq_flag)

    def all_have_flag(self, flag):
        return all(flag in one for one in self._seq_flag)

    def pos2index(self, pos, first):
        if first == 'first':
            return self._seq_pos.index(pos)
        elif first == 'last':
            return len(self._seq_pos) - self._seq_pos[-1::-1].index(pos) - 1
        else:
            raise ValueError('invalid second parameter')

    def posrange2indexrange(self, pos_start, pos_end):
        if self.empty:
            return 0, 0
        if pos_start > self.max_pos:
            idx_start = idx_end = self.pos2index(self.max_pos, 'last')
        elif pos_end < self.min_pos:
            idx_start = idx_end = self.pos2index(self.min_pos, 'first')
        else:
            pos_start = max(self.min_pos, pos_start)
            pos_end = min(self.max_pos, pos_end)
            for pos in range(pos_start, pos_end + 1):
                try:
                    idx_start = self.pos2index(pos, 'first')
                    break
                except ValueError:
                    continue
            for pos in range(pos_end, pos_start - 1, -1):
                try:
                    idx_end = self.pos2index(pos, 'last') + 1
                    break
                except ValueError:
                    continue
        return idx_start, idx_end

    def __getitem__(self, index):
        return PositionalSeqStr(
            self._seq_text[index],
            self._seq_pos[index],
            self._seq_flag[index])

    def __setitem__(self, index, value):
        if isinstance(value, PositionalSeqStr):
            self._seq_text[index] = value._seq_text
            self._seq_pos[index] = value._seq_pos
            self._seq_flag[index] = value._seq_flag
        elif not value:
            self._seq_text[index] = ''
            self._seq_pos[index] = []
            self._seq_flag[index] = []
        elif all(na in GAP_CHARS for na in value):
            self._seq_text[index] = value
            self._seq_pos[index] = [-1] * len(value)
            self._seq_flag[index] = [set() for _ in value]
        else:
            raise ValueError('invalid input')

    def __len__(self):
        return len(self._seq_text)

    def __contains__(self, value):
        return value in self._seq_text

    def __str__(self):
        return self._seq_text

    def __repr__(self):
        return self._seq_text

    def __iter__(self):
        for na, pos, flag in zip(self._seq_text,
                                 self._seq_pos,
                                 self._seq_flag):
            yield PositionalSeqStr(na, pos, flag)

    def __iadd__(self, other):
        self._seq_text += other._seq_text
        self._seq_pos += other._seq_pos
        self._seq_flag += other._seq_flag

    def __add__(self, other):
        seq_text = self._seq_text + other._seq_text
        seq_pos = self._seq_pos + other._seq_pos
        seq_flag = self._seq_flag + other._seq_flag
        return PositionalSeqStr(seq_text, seq_pos, seq_flag)

    def count(self, *args, **kwargs):
        return self._seq_text.count(*args, **kwargs)

    def remove_gaps(self):
        return PositionalSeqStr([
            (na, pos, flag)
            for na, pos, flag in zip(self._seq_text,
                                     self._seq_pos,
                                     self._seq_flag)
            if na not in GAP_CHARS
        ])

    def is_gap(self):
        return bool(
            self._seq_text and
            all(na in GAP_CHARS for na in self._seq_text)
        )

    @classmethod
    def join(cls, value):
        seq_text = ''.join(one._seq_text for one in value)
        seq_pos = [one._seq_pos for one in value]
        return cls(seq_text, seq_pos)


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
        if skip_invalid is not SKIP_VALIDATION:
            seqtext = PositionalSeqStr(seqtext)
            testseqtext = str(seqtext).upper()
            illegal_pattern = ILLEGAL_PATTERNS[seqtype]
            invalids = illegal_pattern.findall(testseqtext)
            if invalids and skip_invalid:
                seqtext = illegal_pattern.sub('-', testseqtext)
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
    def headerdesc(self):
        hd = self.header
        if self.description:
            hd += ' ' + self.description
        return hd

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
            description=self.description,
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
