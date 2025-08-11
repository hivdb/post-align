"""Tests for PAF parser helpers."""

from __future__ import annotations


def test_insert_unaligned_region_negative_seqsize_noop() -> None:
    """Negative unaligned seq size should leave texts unchanged."""
    from postalign.parsers.paf import insert_unaligned_region
    from postalign.models import NAPosition

    reftext = NAPosition.init_from_bytes(b"AAAA")
    seqtext = NAPosition.init_from_bytes(b"AAAA")
    orig = seqtext[:]
    insert_unaligned_region(
        reftext,
        seqtext,
        orig,
        NAPosition,
        align1_ref_end=2,
        align1_seq_end=2,
        align2_ref_start=3,
        align2_seq_start=1,
    )
    assert NAPosition.as_str(reftext) == "AAAA"
    assert NAPosition.as_str(seqtext) == "AAAA"


def test_insert_unaligned_region_inserts_near_alignment2() -> None:
    """Insertion close to alignment2 should offset toward the end."""
    from postalign.parsers.paf import insert_unaligned_region
    from postalign.models import NAPosition

    reftext = NAPosition.init_from_bytes(b"ABCD")
    seqtext = NAPosition.init_from_bytes(b"ABCD")
    orig = NAPosition.init_from_bytes(b"ABXXCD")
    insert_unaligned_region(
        reftext,
        seqtext,
        orig,
        NAPosition,
        align1_ref_end=2,
        align1_seq_end=2,
        align2_ref_start=4,
        align2_seq_start=4,
        insert_close_to=2,
    )
    assert NAPosition.as_str(reftext) == "AB--CD"
    assert NAPosition.as_str(seqtext) == "ABXXCD"
