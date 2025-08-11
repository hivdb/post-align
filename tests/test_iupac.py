"""Tests for IUPAC utilities."""

from __future__ import annotations

from typing import Callable, cast

import pytest


@pytest.fixture()
def iupac_score_func() -> Callable[[int, int], float]:
    """Return the :func:`~postalign.utils.iupac.iupac_score` function."""
    from postalign.utils.iupac import iupac_score

    return cast(Callable[[int, int], float], iupac_score)


def test_iupac_score_match(
    iupac_score_func: Callable[[int, int], float]
) -> None:
    """Identical nucleotides should score 1."""
    assert iupac_score_func(ord(b"A"), ord(b"A")) == 1


def test_iupac_score_deletion(
    iupac_score_func: Callable[[int, int], float]
) -> None:
    """In-frame deletions should score 0."""
    gap = ord(b"-")
    assert iupac_score_func(gap, gap) == 0


def test_iupac_score_ambiguous(
    iupac_score_func: Callable[[int, int], float]
) -> None:
    """Ambiguous mismatch should yield negative fractional score."""
    score = iupac_score_func(ord(b"W"), ord(b"A"))
    assert score == -0.5
