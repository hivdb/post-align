"""Tests for the package version module."""

from __future__ import annotations

def test_version_string() -> None:
    """The version constant should be a non-empty string."""
    from postalign import version

    assert isinstance(version.VERSION, str)
    assert version.VERSION
