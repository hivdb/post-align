"""Tests for the :mod:`postalign.entry` module."""

from __future__ import annotations


def test_entry_exports_cli_and_processors() -> None:
    """Entry module should expose CLI app and processors package."""
    import postalign.entry as entry
    import postalign.cli as cli_module
    import postalign.processors as processors_module

    assert entry.cli is cli_module.cli
    assert entry.processors is processors_module
