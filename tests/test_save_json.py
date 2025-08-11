"""Tests for JSON saving utilities."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest


def test_gene_range_tuples_callback_parses_ranges() -> None:
    """Callback should group gene ranges into tuples."""
    from postalign.processors.save_json import gene_range_tuples_callback

    ctx = MagicMock()
    param = MagicMock()
    result = gene_range_tuples_callback(ctx, param, ("G", "1", "4"))
    assert result == [("G", [(1, 4)])]


def test_gene_range_tuples_callback_missing_value() -> None:
    """Odd range values should raise :class:`BadParameter`."""
    from postalign.processors.save_json import gene_range_tuples_callback
    import typer

    ctx = MagicMock()
    param = MagicMock()
    with pytest.raises(typer.BadParameter):
        gene_range_tuples_callback(ctx, param, ("G", "1"))


def test_gene_range_tuples_callback_refstart_too_small() -> None:
    """Reference start below 1 should raise :class:`BadParameter`."""
    from postalign.processors.save_json import gene_range_tuples_callback
    import typer

    ctx = MagicMock()
    param = MagicMock()
    with pytest.raises(typer.BadParameter):
        gene_range_tuples_callback(ctx, param, ("G", "0", "4"))


def test_gene_range_tuples_callback_refend_too_small() -> None:
    """Reference end must allow at least one codon."""
    from postalign.processors.save_json import gene_range_tuples_callback
    import typer

    ctx = MagicMock()
    param = MagicMock()
    with pytest.raises(typer.BadParameter):
        gene_range_tuples_callback(ctx, param, ("G", "2", "3"))
