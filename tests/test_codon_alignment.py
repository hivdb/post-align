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


def test_extend_codons_until_gap_right() -> None:
    """Scanning forward should stop before codons containing gaps."""
    from postalign.models import NAPosition
    from postalign.processors.codon_alignment import (
        extend_codons_until_gap,
        RIGHT,
    )

    codon1 = NAPosition.init_from_bytes(b"ATG")
    codon2 = NAPosition.init_from_bytes(b"A-G")
    ref_cd, seq_cd, length = extend_codons_until_gap(
        [codon1, codon2], [codon1, codon2], RIGHT
    )
    assert length == 1
    assert [NAPosition.as_str(c) for c in ref_cd] == ["ATG"]


def test_find_windows_with_gap() -> None:
    """Gap windows should be separated when distance exceeds threshold."""
    from postalign.models import NAPosition
    from postalign.processors.codon_alignment import find_windows_with_gap

    ref = NAPosition.init_from_bytes(b"AAA---BBB---CCC")
    seq = NAPosition.init_from_bytes(b"AAA---BBB---CCC")
    windows = find_windows_with_gap(ref, seq, 1)
    assert windows == [slice(3, 6), slice(9, 12)]


def test_find_first_gap_none() -> None:
    """Sequences without gaps should return ``-1`` for first gap."""
    from postalign.models import NAPosition
    from postalign.processors.codon_alignment import find_first_gap

    nas = NAPosition.init_from_bytes(b"ATG")
    assert find_first_gap(nas) == -1


def test_move_gap_to_codon_end() -> None:
    """Gaps within codons should be shifted to the end."""
    from postalign.models import NAPosition
    from postalign.processors.codon_alignment import move_gap_to_codon_end

    codons = [NAPosition.init_from_bytes(b"A-A")]
    moved = move_gap_to_codon_end(codons)
    assert [NAPosition.as_str(c) for c in moved] == ["AA-"]


def test_remove_redundant_gaps() -> None:
    """Matched gaps across sequences should be pruned."""
    from postalign.models import NAPosition
    from postalign.processors.codon_alignment import remove_redundant_gaps

    ref = NAPosition.init_from_bytes(b"A--B")
    seq = NAPosition.init_from_bytes(b"A--B")
    new_ref, new_seq = remove_redundant_gaps(ref, seq)
    assert NAPosition.as_str(new_ref) == "AB"
    assert NAPosition.as_str(new_seq) == "AB"
