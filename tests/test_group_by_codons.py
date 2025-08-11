"""Tests for grouping nucleotides by codons and genes."""

from __future__ import annotations


def test_group_by_codons_basic() -> None:
    """Sequences should be split into codons accounting for gaps."""
    from postalign.models import NAPosition
    from postalign.utils.group_by_codons import group_by_codons

    ref = NAPosition.init_from_bytes(b"ATGCTA")
    seq = NAPosition.init_from_bytes(b"ATG-TA")
    ref_codons, seq_codons = group_by_codons(ref, seq)
    assert [NAPosition.as_str(c) for c in ref_codons] == ["ATG", "CTA"]
    assert [NAPosition.as_str(c) for c in seq_codons] == ["ATG", "-TA"]


def test_find_codon_trim_slice() -> None:
    """Leading and trailing gap-only codons should be trimmed."""
    from postalign.models import NAPosition
    from postalign.utils.group_by_codons import find_codon_trim_slice

    codons = [
        NAPosition.init_from_bytes(b"---"),
        NAPosition.init_from_bytes(b"ATG"),
        NAPosition.init_from_bytes(b"---"),
    ]
    trim_slice = find_codon_trim_slice(codons)
    assert trim_slice == slice(1, 2)


def test_group_by_gene_codons() -> None:
    """Codons should be collected per gene range."""
    from postalign.models import NAPosition
    from postalign.utils.group_by_codons import group_by_gene_codons

    ref = NAPosition.init_from_bytes(b"ATGAAA")
    seq = NAPosition.init_from_bytes(b"ATGAAA")
    genes = [("geneA", [(1, 3)]), ("geneB", [(4, 6)])]
    results = group_by_gene_codons(ref, seq, genes)
    assert results[0][0] == "geneA"
    assert [NAPosition.as_str(c) for c in results[0][1]] == ["ATG"]
    assert results[1][0] == "geneB"
    assert [NAPosition.as_str(c) for c in results[1][1]] == ["AAA"]
