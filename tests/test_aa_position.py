"""Tests for the amino acid position placeholder."""

from __future__ import annotations

import pytest


def test_aaposition_not_implemented() -> None:
    """AAPosition methods should raise ``NotImplementedError``."""
    from postalign.models import AAPosition

    with pytest.raises(NotImplementedError):
        AAPosition.init_gaps(1)
