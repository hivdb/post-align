"""Tests for codon alignment utilities."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest
import typer

from postalign.processors.codon_alignment import (
    REFGAP,
    SEQGAP,
    codon_alignment,
    gap_placement_score_callback,
    parse_gap_placement_score,
)


def test_parse_gap_placement_score_valid() -> None:
    """Valid gap placement strings should map to score dictionaries."""

    scores = parse_gap_placement_score("204ins:-5,2041/12del:10")
    assert scores[REFGAP][(204, 0)] == -5
    assert scores[SEQGAP][(2041, 12)] == 10


def test_parse_gap_placement_score_invalid() -> None:
    """Invalid score strings should raise :class:`ValueError`."""

    with pytest.raises(ValueError):
        parse_gap_placement_score("204foo")


def test_gap_placement_score_callback_missing_name() -> None:
    """Missing parameter name should raise :class:`BadParameter`."""

    ctx = MagicMock()
    param = MagicMock()
    param.name = None
    with pytest.raises(typer.BadParameter):
        gap_placement_score_callback(ctx, param, ())


def test_codon_alignment_invalid_ref_start() -> None:
    """Invalid reference start should raise :class:`BadParameter`."""

    with pytest.raises(typer.BadParameter):
        codon_alignment(ref_start=0, ref_end=-1)


def test_codon_alignment_invalid_ref_end() -> None:
    """Too-small ref_end should raise :class:`BadParameter`."""

    with pytest.raises(typer.BadParameter):
        codon_alignment(ref_start=5, ref_end=6)
