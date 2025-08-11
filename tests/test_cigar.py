"""Tests for CIGAR utilities."""

from __future__ import annotations

import pytest
from typing import Any, Type


@pytest.fixture()
def cigar_cls() -> Type[Any]:
    """Return the :class:`~postalign.utils.cigar.CIGAR` class."""
    from postalign.utils.cigar import CIGAR

    return CIGAR


@pytest.fixture()
def na_position_cls() -> Type[Any]:
    """Return the :class:`~postalign.models.NAPosition` class."""
    from postalign.models import NAPosition

    return NAPosition


def test_shrink_by_ref_with_insertion(cigar_cls: Type[Any]) -> None:
    """Insertions should remain when shrinking by reference length."""
    cigar = cigar_cls(0, 0, "5M2I5M")
    shrunk = cigar.shrink_by_ref(6)
    assert shrunk.get_cigar_string() == "5M2I1M"


def test_get_alignment_handles_insertion(
    cigar_cls: Type[Any], na_position_cls: Type[Any]
) -> None:
    """``get_alignment`` should pad reference with gaps for insertions."""
    ref = na_position_cls.init_from_bytes(b"ACGT")
    seq = na_position_cls.init_from_bytes(b"ATCGT")
    cigar = cigar_cls(0, 0, "1M1I3M")
    ref_aln, seq_aln = cigar.get_alignment(ref, seq, na_position_cls)
    assert "".join(str(p) for p in ref_aln) == "A-CGT"
    assert "".join(str(p) for p in seq_aln) == "ATCGT"


def test_get_alignment_handles_deletion(
    cigar_cls: Type[Any], na_position_cls: Type[Any]
) -> None:
    """Deletion operations should insert gaps into the sequence."""
    ref = na_position_cls.init_from_bytes(b"ACGT")
    seq = na_position_cls.init_from_bytes(b"AGT")
    cigar = cigar_cls(0, 0, "1M1D2M")
    ref_aln, seq_aln = cigar.get_alignment(ref, seq, na_position_cls)
    assert "".join(str(p) for p in ref_aln) == "ACGT"
    assert "".join(str(p) for p in seq_aln) == "A-GT"


def test_get_alignment_raises_on_mismatch(
    cigar_cls: Type[Any], na_position_cls: Type[Any]
) -> None:
    """Mismatched alignment lengths should raise ``ValueError``."""
    ref = na_position_cls.init_from_bytes(b"A")
    seq = na_position_cls.init_from_bytes(b"AGC")
    cigar = cigar_cls(0, 0, "1M1I1M")
    with pytest.raises(ValueError):
        cigar.get_alignment(ref, seq, na_position_cls)


def test_cigar_repr(cigar_cls: Type[Any]) -> None:
    """``repr`` should show CIGAR string and offsets."""
    cigar = cigar_cls(1, 2, "5M")
    assert repr(cigar) == "<CIGAR '5M' ref_start=1 seq_start=2>"
