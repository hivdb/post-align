"""Tests for :mod:`postalign.models.na_position`."""

from __future__ import annotations


def test_enumerate_seq_pos_handles_gaps() -> None:
    """Gap characters should be reported with ``-1`` positions."""
    from postalign.models.na_position import enumerate_seq_pos

    assert enumerate_seq_pos(b"A-C") == [1, -1, 2]


def test_min_max_pos_and_init_gaps() -> None:
    """Verify min/max position helpers and gap initialization."""
    from postalign.models.na_position import NAPosition

    nas = NAPosition.init_from_bytes(b"-AC-")
    assert NAPosition.min_pos(nas) == 1
    assert NAPosition.max_pos(nas) == 2
    gaps = NAPosition.init_gaps(2)
    assert all(na.is_gap for na in gaps)
