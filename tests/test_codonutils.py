"""Tests for codon translation utilities."""

from __future__ import annotations


def test_translate_codon_basic() -> None:
    """Standard codon should translate to a single amino acid."""
    from postalign.models import NAPosition
    from postalign.utils.codonutils import translate_codon

    nas = NAPosition.init_from_bytes(b"ATG")
    assert translate_codon(nas) == b"M"


def test_translate_codon_with_ambiguity() -> None:
    """Ambiguous bases should expand to multiple amino acids."""
    from postalign.models import NAPosition
    from postalign.utils.codonutils import translate_codon

    nas = NAPosition.init_from_bytes(b"TGN")
    assert translate_codon(nas) == b"*CW"


def test_translate_codon_inframe_deletion() -> None:
    """A gap-only codon should translate to the deletion marker."""
    from postalign.models import NAPosition
    from postalign.utils.codonutils import translate_codon

    nas = NAPosition.init_gaps(3)
    assert translate_codon(nas, del_as=b"-") == b"-"


def test_translate_codon_frameshift_short() -> None:
    """Codons shorter than three bases should yield frameshift marker."""
    from postalign.models import NAPosition
    from postalign.utils.codonutils import translate_codon

    nas = NAPosition.init_from_bytes(b"AC")
    assert translate_codon(nas, fs_as=b"X") == b"X"


def test_translate_codons_calls_internal() -> None:
    """``translate_codons`` should call ``_translate_codon`` for each chunk."""
    from unittest.mock import patch
    from postalign.models import NAPosition
    from postalign.utils import codonutils

    nas = NAPosition.init_from_bytes(b"ATGATG")
    with patch.object(
        codonutils, "_translate_codon", return_value=(ord("A"),)
    ) as mock_trans:
        result = codonutils.translate_codons(nas)
        assert result == [b"A", b"A"]
        assert mock_trans.call_count == 2


def test_get_codons_returns_expected() -> None:
    """Reverse lookup should list codons for a given amino acid."""
    from postalign.utils.codonutils import get_codons

    assert get_codons(ord(b"M")) == [b"ATG"]


def test_compare_codon_match_and_mismatch() -> None:
    """Compare codons allowing for ambiguity."""
    from postalign.utils.codonutils import compare_codon

    assert compare_codon(b"ATG", b"ATG") is True
    assert compare_codon(b"ATG", b"ATR") is True
    assert compare_codon(b"ATG", b"ATN") is False
    assert compare_codon(b"ATG", b"ATA") is False


def test_compare_codon_ambiguous_mismatch() -> None:
    """Mismatched ambiguous bases should return ``False``."""
    from postalign.utils.codonutils import compare_codon

    assert compare_codon(b"ACG", b"ACY") is False
