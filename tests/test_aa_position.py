"""Tests for the amino acid position placeholder."""

from __future__ import annotations

import pytest


def test_aaposition_not_implemented() -> None:
    """AAPosition methods should raise ``NotImplementedError``."""
    from postalign.models import AAPosition

    with pytest.raises(NotImplementedError):
        AAPosition.init_gaps(1)


def test_aaposition_any_has_gap_not_implemented() -> None:
    """Gap queries should be unimplemented for AAPosition."""
    from postalign.models import AAPosition

    with pytest.raises(NotImplementedError):
        AAPosition.any_has_gap([])


def test_aaposition_set_flag_not_implemented() -> None:
    """Flag setters should raise :class:`NotImplementedError`."""
    from postalign.models import AAPosition
    from postalign.models.position_flag import PositionFlag

    with pytest.raises(NotImplementedError):
        AAPosition.set_flag([], PositionFlag.NONE)
