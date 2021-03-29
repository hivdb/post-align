import re
from collections import namedtuple

from .modifier import ModifierLinkedList

GAP_CHARS = '.-'
GAP_PATTERN = re.compile(r'^[.-]+$')
ANY_GAP_PATTERN = re.compile(r'[.-]')

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

    @classmethod
    def init_empty(cls):
        return cls('', [], [])

    @classmethod
    def init_gaps(cls, gaplen):
        seq_text = '-' * gaplen
        seq_pos = [-1] * gaplen
        seq_flag = [set() for _ in seq_pos]
        return cls(seq_text, seq_pos, seq_flag)

    @classmethod
    def init_from_triplets(cls, triplets):
        (seq_text,
         seq_pos,
         seq_flag) = list(zip(*triplets))
        return cls(seq_text, seq_pos, seq_flag)

    @classmethod
    def init_from_nastring(cls, seq_text):
        seq_text = seq_text.upper()
        seq_pos = enumerate_seq_pos(seq_text)
        seq_flag = [set() for _ in seq_pos]
        return cls(seq_text, seq_pos, seq_flag)

    def __init__(self, seq_text, seq_pos, seq_flag):
        self._seq_text = seq_text
        self._seq_pos = seq_pos
        self._seq_flag = seq_flag
        self._min_pos = None
        self._max_pos = None
        if len(seq_text) == 1:
            self._is_gap = seq_text in GAP_CHARS
            self.is_single_gap = self._is_gap
        else:
            self._is_gap = None
            self.is_single_gap = False

    @property
    def min_pos(self):
        if self._min_pos is None:
            self._min_pos = next(
                (pos for pos in self._seq_pos if pos > 0), None
            )
        return self._min_pos

    @property
    def max_pos(self):
        if self._max_pos is None:
            self._max_pos = next(
                (pos for pos in reversed(self._seq_pos) if pos > 0), None
            )
        return self._max_pos

    @property
    def empty(self):
        return self.min_pos is None

    def set_flag(self, flag):
        for one in self._seq_flag:
            one.add(flag)

    def any_has_flag(self, flag):
        for one in self._seq_flag:
            if flag in one:
                return True
        return False

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
        if isinstance(index, slice):
            return PositionalSeqStr(
                self._seq_text[index],
                self._seq_pos[index],
                self._seq_flag[index])
        else:
            return PositionalSeqStr(
                self._seq_text[index],
                [self._seq_pos[index]],
                [self._seq_flag[index]])

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
            yield PositionalSeqStr(na, [pos], [flag])

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
        return PositionalSeqStr.init_from_triplets([
            (na, pos, flag)
            for na, pos, flag in zip(self._seq_text,
                                     self._seq_pos,
                                     self._seq_flag)
            if na not in GAP_CHARS
        ])

    def is_gap(self):
        if self._is_gap is None:
            self._is_gap = bool(
                self._seq_text and
                GAP_PATTERN.match(self._seq_text)
            )
        return self._is_gap

    @classmethod
    def join(cls, value):
        seq_text = ''.join(one._seq_text for one in value)
        seq_pos = [pos for one in value for pos in one._seq_pos]
        seq_flag = [flag for one in value for flag in one._seq_flag]
        return cls(seq_text, seq_pos, seq_flag)

    @staticmethod
    def list_contains_any_gap(nas):
        for na in nas:
            if na.is_gap():
                return True
        return False

    @staticmethod
    def list_contains_all_gap(nas):
        for na in nas:
            if not na.is_gap():
                return False
        return True


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
        if not isinstance(seqtext, PositionalSeqStr):
            raise ValueError(
                'seqtext must be an instance of PositionSeqStr'
            )
        if skip_invalid is not SKIP_VALIDATION:
            testseqtext = str(seqtext)
            illegal_pattern = ILLEGAL_PATTERNS[seqtype]
            invalids = illegal_pattern.findall(testseqtext)
            if invalids and skip_invalid:
                seqtext = PositionalSeqStr.init_from_nastring(
                    illegal_pattern.sub('-', testseqtext)
                )
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
