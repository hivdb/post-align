import re
import weakref
from itertools import groupby
from collections import namedtuple


SEQTYPES = ['NA']  # currently we don't support AA sequence

SKIP_VALIDATION = object()

ILLEGAL_PATTERNS = {
    'NA': re.compile(r'[^ACGTUWSMKRYBDHVN.-]')
}

BaseSequence = namedtuple(
    'BaseSequence', ['header', 'seqtext', 'seqid', 'seqtype', 'modifiers_']
)


class Modifier:

    def __init__(self, text, *, slicetuples=None):
        self.text = text
        self.slicetuples = slicetuples or []
        self.child_mods = set()
        self.parent_mods = []
        is_root = text == 'root()'
        self.step = 0 if is_root else None
        self.root_modifier = self

    def add_child_mod(self, mod):
        self.child_mods.add(mod)
        mod.parent_mods.append(weakref.ref(self))
        mod.root_modifier = self.root_modifier
        mod.step = self.step + 1

    def remove_child_mod(self, mod):
        self.child_mods.remove(mod)

    def get_all_offspring_mods(self):
        offspring_mods = set()
        for mod in self.child_mods:
            offspring_mods.add(mod)
            offspring_mods |= mod.get_all_offspring_mods()
        return offspring_mods

    def __add__(self, other):
        if not isinstance(other, Modifier):
            raise TypeError(
                'unsupported operand type(s) for +: {!r} and {!r}'
                .format(self.__class__.__name__, other.__class__.__name__)
            )
        if (
            self.slicetuples and other.slicetuples and
            len(self.parent_mods) == len(other.parent_mods) and
            self.parent_mods[0]() == other.parent_mods[0]()
        ):
            # special case compress modifier: two slice + join
            slicetuples = self.slicetuples + other.slicetuples
            modtext = 'join({})'.format(
                ','.join('{}..{}'.format(*st) for st in slicetuples)
            )
            merged = Modifier(modtext, slicetuples=slicetuples)
            parent = self.parent_mods[0]()
            parent.remove_child_mod(self)
            parent.remove_child_mod(other)
            parent.add_child_mod(merged)
        else:
            merged = Modifier('concat(step{},step{})'
                              .format(self.step, other.step))
            self.add_child_mod(merged)
            other.add_child_mod(merged)
            merged.root_modifier = self.root_modifier
            merged.step = max(self.step, other.step) + 1
        return merged

    def __str__(self):
        return self.text


class ModifierLinkedList:

    def __init__(self, last_modifier=None):
        if last_modifier is None:
            last_modifier = Modifier('root()')
        self._last_modifier = last_modifier

    def push(self, modtext, **kw):
        modifier = Modifier(modtext, **kw)
        self._last_modifier.add_child_mod(modifier)
        return ModifierLinkedList(modifier)

    def replace_last(self, modtext, **kw):
        modifier = Modifier(modtext, **kw)
        last_modifier = self._last_modifier
        for parent in last_modifier.parent_mods:
            parent = parent()
            parent.remove_child_mod(last_modifier)
            parent.add_child_mod(modifier)
        return ModifierLinkedList(modifier)

    @property
    def last_modifier(self):
        return self._last_modifier

    def __add__(self, other):
        if not isinstance(other, ModifierLinkedList):
            raise TypeError(
                'unsupported operand type(s) for +: {!r} and {!r}'
                .format(self.__class__.__name__, other.__class__.__name__)
            )
        last_modifier = self._last_modifier + other._last_modifier
        return ModifierLinkedList(last_modifier)

    def __str__(self):
        root = self.last_modifier.root_modifier
        all_mods = root.get_all_offspring_mods()
        grouped_mods = groupby(
            sorted(all_mods, key=lambda mod: mod.step),
            lambda mod: mod.step
        )
        text_buffer = []
        for step, group in grouped_mods:
            text_buffer.append('{}:{}'.format(
                step, ','.join(str(m) for m in group)))
        return '|'.join(text_buffer)


class Sequence(BaseSequence):

    def __new__(
        cls, *, header, seqtext, seqid,
        seqtype, modifiers_=None, skip_invalid=True
    ):
        illegal_pattern = ILLEGAL_PATTERNS[seqtype]
        invalids = illegal_pattern.findall(seqtext)
        if skip_invalid is not SKIP_VALIDATION:
            if invalids and skip_invalid:
                seqtext = illegal_pattern.sub('.', seqtext)
            elif invalids:
                raise ValueError(
                    'sequence {} contains invalid notation(s) ({})'
                    'while skip_invalid=False'
                    .format(header, ''.join(set(invalids)))
                )
        return super().__new__(cls, header, seqtext,
                               seqid, seqtype, modifiers_)

    @property
    def modifiers(self):
        if self.modifiers_ is None:
            return ModifierLinkedList()
        else:
            return self.modifiers_

    def __getitem__(self, index):
        seqtext = self.seqtext[index]
        if isinstance(index, slice):
            sliceval = index.indices(len(self.seqtext))
            if sliceval[2] == 1:
                modtext = 'slice({},{})'.format(*sliceval[:2])
            else:
                raise ValueError(
                    'step slicing is unsupported: {!r}'.format(index)
                )
            modifiers = self.modifiers.push(
                modtext, slicetuples=[sliceval[:2]])

            return Sequence(
                header=self.header,
                seqtext=seqtext,
                seqid=self.seqid,
                seqtype=self.seqtype,
                modifiers_=modifiers,
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
            skip_invalid=SKIP_VALIDATION)

    def __iter__(self):
        yield from self.seqtext

    def __len__(self):
        return len(self.seqtext)

    def push_seqtext(self, seqtext, modtext, **kw):
        """Modify seqtext and push modifier forward

        The difference between `push_seqtext` and `replace_seqtext` is similar
        to modern browsers' `history.pushstate` and `history.replacestate`,
        where the history is stored in attribute `modifiers`.
        """
        modifiers = self.modifiers.push(modtext, **kw)

        return Sequence(
            header=self.header,
            seqtext=seqtext,
            seqid=self.seqid,
            seqtype=self.seqtype,
            modifiers_=modifiers,
            skip_invalid=True)

    def replace_seqtext(self, seqtext, modtext, **kw):
        """Modify seqtext and replace last modifier

        The difference between `push_seqtext` and `replace_seqtext` is similar
        to modern browsers' `history.pushstate` and `history.replacestate`,
        where the history is stored in attribute `modifiers`.
        """
        modifiers = self.modifiers.replace_last(modtext, **kw)

        return Sequence(
            header=self.header,
            seqtext=seqtext,
            seqid=self.seqid,
            seqtype=self.seqtype,
            modifiers_=modifiers,
            skip_invalid=True)

    @property
    def header_with_modifiers(self):
        return '{} MOD::{}'.format(self.header, self.modifiers)
