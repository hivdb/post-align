"""Tests for sequence sanitization helpers."""

from __future__ import annotations

import pytest


def test_sanitize_sequence_skips_invalid() -> None:
    """Invalid nucleotides should be removed when skipping is enabled."""
    from postalign.models.na_position import NAPosition
    from postalign.models._sequence import sanitize_sequence

    nas = NAPosition.init_from_bytes(b"AZT")
    cleaned = sanitize_sequence(nas, NAPosition, "hdr", True)
    assert [str(na) for na in cleaned] == ["A", "T"]


def test_sanitize_sequence_raises_on_invalid() -> None:
    """Disallowed nucleotides should raise an error when not skipped."""
    from postalign.models.na_position import NAPosition
    from postalign.models._sequence import sanitize_sequence

    nas = NAPosition.init_from_bytes(b"AZT")
    with pytest.raises(ValueError):
        sanitize_sequence(nas, NAPosition, "hdr", False)
