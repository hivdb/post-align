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


def test_translate_codon_stop() -> None:
    """Stop codons should translate to the stop symbol."""
    from postalign.models import NAPosition
    from postalign.utils.codonutils import translate_codon

    nas = NAPosition.init_from_bytes(b"TGA")
    assert translate_codon(nas) == b"*"


def test_translate_codons_translates_chunks() -> None:
    """``translate_codons`` should translate each codon chunk."""
    from unittest.mock import patch
    from postalign.models import NAPosition
    from postalign.utils import codonutils

    nas = NAPosition.init_from_bytes(b"ATGTTT")

    from collections.abc import Iterator

    def chunk_list(seq: list, n: int) -> Iterator[list]:
        for i in range(0, len(seq), n):
            yield list(seq[i:i + n])

    with patch.object(codonutils, "chunked", chunk_list):
        assert codonutils.translate_codons(nas) == [b"M", b"F"]


def test_get_codons_returns_expected() -> None:
    """Reverse lookup should list codons for a given amino acid."""
    from postalign.utils.codonutils import get_codons

    assert get_codons(ord(b"M")) == [b"ATG"]


def test_get_codons_unknown_raises_key_error() -> None:
    """Unknown amino acids should raise :class:`KeyError`."""
    from postalign.utils.codonutils import get_codons
    import pytest

    with pytest.raises(KeyError):
        get_codons(ord(b"?"))


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
