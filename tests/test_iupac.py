"""Tests for IUPAC utilities."""

from postalign.utils.iupac import iupac_score


def test_iupac_score_match() -> None:
    """Identical nucleotides should score 1."""
    assert iupac_score(ord(b'A'), ord(b'A')) == 1


def test_iupac_score_deletion() -> None:
    """In-frame deletions should score 0."""
    gap = ord(b'-')
    assert iupac_score(gap, gap) == 0


def test_iupac_score_ambiguous() -> None:
    """Ambiguous mismatch should yield negative fractional score."""
    score = iupac_score(ord(b'W'), ord(b'A'))
    assert score == -0.5
