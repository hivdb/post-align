"""Tests for IUPAC nucleotide scoring."""


def test_iupac_score_gap_match() -> None:
    """Matching gaps should incur no penalty."""
    from postalign.utils.iupac import iupac_score

    gap = ord(b"-")
    assert iupac_score(gap, gap) == 0


def test_iupac_score_partial_mismatch() -> None:
    """Ambiguous bases should yield a fractional mismatch score."""
    from postalign.utils.iupac import iupac_score

    assert iupac_score(ord(b"W"), ord(b"A")) == -0.5
