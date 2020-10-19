import weakref
from itertools import groupby
from collections import Counter


class Modifier:

    def __init__(self, text, *, slicetuples=None):
        self.text = text
        self.slicetuples = slicetuples or []
        self.child_mods = Counter()
        self.parent_mods = []
        is_root = text == 'root()'
        self.step = 0 if is_root else None
        self.root_modifier = self

    def add_child_mod(self, mod):
        self.child_mods[mod] += 1
        mod.parent_mods.append(weakref.ref(self))
        mod.root_modifier = self.root_modifier
        mod.step = self.step + 1

    def remove_child_mod(self, mod):
        if mod not in self.child_mods:
            raise KeyError('Modifier {} not found'.format(mod))
        self.child_mods[mod] -= 1

    def get_all_offspring_mods(self):
        offspring_mods = set()
        for mod, count in self.child_mods.items():
            if count <= 0:
                continue  # ignore removed child_mod
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
        text_buffer = []
        for step, group in self:
            text_buffer.append('{}:{}'.format(
                step, ','.join(str(m) for m in group)))
        return '|'.join(text_buffer)

    def __iter__(self):
        root = self.last_modifier.root_modifier
        all_mods = root.get_all_offspring_mods()
        return groupby(
            sorted(all_mods, key=lambda mod: mod.step),
            lambda mod: mod.step
        )
