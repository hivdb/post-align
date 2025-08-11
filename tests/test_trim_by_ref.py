"""Tests for trimming sequences by reference gaps."""

from __future__ import annotations

from postalign.models import Sequence, NAPosition


def _make_seq(header: str, text: bytes, seqid: int) -> Sequence:
    """Construct a sequence with the given text and identifier."""

    seqtext = NAPosition.init_from_bytes(text)
    return Sequence(
        header=header,
        description="",
        seqtext=seqtext,
        seqid=seqid,
        seqtype=NAPosition,
        abs_seqstart=1,
        skip_invalid=True,
    )


def test_find_trim_slice() -> None:
    """Leading/trailing gaps should be reflected in the slice."""
    from postalign.processors.trim_by_ref import find_trim_slice

    ref = _make_seq("ref", b"--ACGT--", 1)
    assert find_trim_slice(ref) == slice(2, 6)


def test_trim_by_ref_processor_trims() -> None:
    """Processor should trim alignments based on reference gaps."""
    from postalign.processors.trim_by_ref import trim_by_ref

    ref = _make_seq("ref", b"--ACGT--", 1)
    seq = _make_seq("seq", b"--AC-T--", 2)
    proc = trim_by_ref()
    trimmed_ref, trimmed_seq = list(proc([(ref, seq)], []))[0]
    assert trimmed_ref.seqtext_as_str == "ACGT"
    assert trimmed_seq.seqtext_as_str == "AC-T"
