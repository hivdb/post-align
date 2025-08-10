"""Tests for CIGAR utilities."""

import pytest

from postalign.utils.cigar import CIGAR
from postalign.models import NAPosition


def test_shrink_by_ref_with_insertion() -> None:
    """Insertions should remain when shrinking by reference length."""
    cigar = CIGAR(0, 0, "5M2I5M")
    shrunk = cigar.shrink_by_ref(6)
    assert shrunk.get_cigar_string() == "5M2I1M"


def test_get_alignment_handles_insertion() -> None:
    """``get_alignment`` should pad reference with gaps for insertions."""
    ref = NAPosition.init_from_bytes(b"ACGT")
    seq = NAPosition.init_from_bytes(b"ATCGT")
    cigar = CIGAR(0, 0, "1M1I3M")
    ref_aln, seq_aln = cigar.get_alignment(ref, seq, NAPosition)
    assert "".join(str(p) for p in ref_aln) == "A-CGT"
    assert "".join(str(p) for p in seq_aln) == "ATCGT"


def test_get_alignment_raises_on_mismatch() -> None:
    """Mismatched alignment lengths should raise ``ValueError``."""
    ref = NAPosition.init_from_bytes(b"A")
    seq = NAPosition.init_from_bytes(b"AGC")
    cigar = CIGAR(0, 0, "1M1I1M")
    with pytest.raises(ValueError):
        cigar.get_alignment(ref, seq, NAPosition)
