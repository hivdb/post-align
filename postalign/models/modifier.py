import weakref
from operator import attrgetter
from typing import (
    Any, List, Set, Tuple, Iterable, Iterator, Optional, Callable
)
from itertools import groupby
from collections import Counter


class Modifier:

    text: str
    slicetuples: List[Tuple[int, int]]
    child_mods: Counter['Modifier']
    parent_mods: List[weakref.ref['Modifier']]
    step: Optional[int]
    root_modifier: 'Modifier'

    def __init__(
        self: 'Modifier',
        text: str,
        *,
        slicetuples: Optional[List[Tuple[int, int]]] = None
    ) -> None:
        self.text = text
        self.slicetuples = slicetuples or []
        self.child_mods = Counter()
        self.parent_mods = []
        is_root: bool = text == 'root()'
        self.step = 0 if is_root else None
        self.root_modifier = self

    def add_child_mod(self: 'Modifier', mod: 'Modifier') -> None:
        self.child_mods[mod] += 1
        mod.parent_mods.append(weakref.ref(self))
        mod.root_modifier = self.root_modifier
        if self.step is not None:
            mod.step = self.step + 1

    def remove_child_mod(self: 'Modifier', mod: 'Modifier') -> None:
        if mod not in self.child_mods:
            raise KeyError('Modifier {} not found'.format(mod))
        self.child_mods[mod] -= 1

    def get_all_offspring_mods(self: 'Modifier') -> Set['Modifier']:
        mod: Modifier
        count: int
        offspring_mods: Set[Modifier] = set()
        for mod, count in self.child_mods.items():
            if count <= 0:
                continue  # ignore removed child_mod
            offspring_mods.add(mod)
            offspring_mods |= mod.get_all_offspring_mods()
        return offspring_mods

    def __add__(self: 'Modifier', other: 'Modifier') -> 'Modifier':
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
            parent: Optional[Modifier] = self.parent_mods[0]()
            if parent:
                parent.remove_child_mod(self)
                parent.remove_child_mod(other)
                parent.add_child_mod(merged)
        else:
            merged = Modifier('concat(step{},step{})'
                              .format(self.step, other.step))
            self.add_child_mod(merged)
            other.add_child_mod(merged)
            merged.root_modifier = self.root_modifier
            if self.step and other.step:
                merged.step = max(self.step, other.step) + 1
        return merged

    def __str__(self: 'Modifier') -> str:
        return self.text


class ModifierLinkedList:

    def __init__(
        self: 'ModifierLinkedList',
        last_modifier: Optional[Modifier] = None
    ) -> None:
        if last_modifier is None:
            last_modifier = Modifier('root()')
        self._last_modifier = last_modifier

    def push(
        self: 'ModifierLinkedList',
        modtext: str,
        **kw: Any
    ) -> 'ModifierLinkedList':
        modifier: Modifier = Modifier(modtext, **kw)
        self._last_modifier.add_child_mod(modifier)
        return ModifierLinkedList(modifier)

    def replace_last(
        self: 'ModifierLinkedList',
        modtext: str,
        **kw: Any
    ) -> 'ModifierLinkedList':
        parent_ref: weakref.ref[Modifier]
        modifier: Modifier = Modifier(modtext, **kw)
        last_modifier: Modifier = self._last_modifier
        for parent_ref in last_modifier.parent_mods:
            parent: Optional[Modifier] = parent_ref()
            if parent:
                parent.remove_child_mod(last_modifier)
                parent.add_child_mod(modifier)
        return ModifierLinkedList(modifier)

    @property
    def last_modifier(self: 'ModifierLinkedList') -> Modifier:
        return self._last_modifier

    def __add__(
        self: 'ModifierLinkedList',
        other: 'ModifierLinkedList'
    ) -> 'ModifierLinkedList':
        if not isinstance(other, ModifierLinkedList):
            raise TypeError(
                'unsupported operand type(s) for +: {!r} and {!r}'
                .format(self.__class__.__name__, other.__class__.__name__)
            )
        last_modifier: Modifier = self._last_modifier + other._last_modifier
        return ModifierLinkedList(last_modifier)

    def __str__(self: 'ModifierLinkedList') -> str:
        step: int
        group: Iterable[Modifier]
        text_buffer: List[str] = []
        for step, group in self:
            text_buffer.append('{}:{}'.format(
                step, ','.join(str(m) for m in group)))
        return '|'.join(text_buffer)

    def __iter__(
        self: 'ModifierLinkedList'
    ) -> Iterator[Tuple[int, Iterable[Modifier]]]:
        root = self.last_modifier.root_modifier
        all_mods = root.get_all_offspring_mods()
        get_step: Callable[[Modifier], int] = attrgetter('step')
        return groupby(
            sorted(all_mods, key=get_step),
            get_step
        )
