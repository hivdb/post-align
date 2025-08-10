"""Tests for the :mod:`postalign.utils.blosum62` module."""

from __future__ import annotations


def test_blosum62_identical() -> None:
    """Identical amino acids should yield their matrix score."""
    from postalign.utils.blosum62 import blosum62_score

    assert blosum62_score(b"A", b"A") == 4


def test_blosum62_deletions_and_frameshift() -> None:
    """Deletion and frameshift markers should apply their penalties."""
    from postalign.utils.blosum62 import blosum62_score

    assert blosum62_score(b"-", b"A") == -1
    assert blosum62_score(b"-", b"-") == 0
    assert blosum62_score(b"X", b"A") == -1


def test_blosum62_unknown_amino_acid() -> None:
    """Unrecognized amino acids contribute zero score."""
    from postalign.utils.blosum62 import blosum62_score

    assert blosum62_score(b"Z", b"Z") == 0


def test_blosum62_empty_sequences() -> None:
    """Empty inputs should yield a zero score without division by zero."""
    from postalign.utils.blosum62 import blosum62_score

    assert blosum62_score(b"", b"") == 0
