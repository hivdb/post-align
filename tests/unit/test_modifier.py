"""Tests for :mod:`postalign.models.modifier`."""

from __future__ import annotations

import pytest


def test_add_and_remove_child_mod() -> None:
    """Adding and removing child modifiers adjusts counts and steps."""
    from postalign.models.modifier import Modifier

    parent = Modifier('root()')
    child = Modifier('child')
    parent.add_child_mod(child)
    assert parent.child_mods[child] == 1
    assert child.step == 1
    parent.remove_child_mod(child)
    assert parent.child_mods[child] == 0


def test_remove_child_mod_missing_raises() -> None:
    """Removing a non-existent child should raise :class:`KeyError`."""
    from postalign.models.modifier import Modifier

    parent = Modifier('root()')
    with pytest.raises(KeyError):
        parent.remove_child_mod(Modifier('missing'))


def test_modifier_merge_slice_join() -> None:
    """Merging sibling slice modifiers should produce a join modifier."""
    from postalign.models.modifier import Modifier

    parent = Modifier('root()')
    mod1 = Modifier('slice1', slicetuples=[(0, 1)])
    mod2 = Modifier('slice2', slicetuples=[(1, 2)])
    parent.add_child_mod(mod1)
    parent.add_child_mod(mod2)
    merged = mod1 + mod2
    assert merged.text.startswith('join(')
    assert parent.child_mods[mod1] == 0
    assert parent.child_mods[mod2] == 0
    assert parent.child_mods[merged] == 1


def test_modifier_add_non_modifier_type_error() -> None:
    """Adding a non-Modifier object should raise :class:`TypeError`."""
    from postalign.models.modifier import Modifier

    mod = Modifier('foo')
    with pytest.raises(TypeError):
        _ = mod + 1  # type: ignore[operator]


def test_modifier_linked_list_push_and_str() -> None:
    """Pushing to the linked list should update its string representation."""
    from postalign.models.modifier import ModifierLinkedList

    root_list = ModifierLinkedList()
    ll = root_list.push('m1')
    assert 'm1' in str(ll)
    step, mods = next(iter(ll))
    assert step == 1
    assert [m.text for m in mods] == ['m1']


def test_modifier_linked_list_replace_last() -> None:
    """Replacing the last modifier should update parent links."""
    from postalign.models.modifier import ModifierLinkedList

    root_list = ModifierLinkedList()
    ll = root_list.push('m1')
    ll2 = ll.replace_last('m2')
    assert ll2.last_modifier.text == 'm2'
    parent = ll2.last_modifier.parent_mods[0]()
    assert parent is not None and parent.child_mods[ll2.last_modifier] == 1


def test_modifier_linked_list_add_lists() -> None:
    """Adding two lists should join their last modifiers."""
    from postalign.models.modifier import ModifierLinkedList

    left = ModifierLinkedList().push('left')
    right = ModifierLinkedList().push('right')
    merged = left + right
    assert merged.last_modifier.text.startswith('concat(')
