"""Tests for the save_fasta processor."""

from __future__ import annotations

def _make_seq(header: str, seqid: int, text: bytes):
    """Build a :class:`Sequence` from raw bytes."""
    from postalign.models import Sequence, NAPosition

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


def test_save_fasta_pairwise_outputs_ref_each_time() -> None:
    """Pairwise mode should output references for every pair."""
    from postalign.processors.save_fasta import save_fasta

    ref1 = _make_seq("ref1", 1, b"AA")
    seq1 = _make_seq("seq1", 2, b"AA")
    ref2 = _make_seq("ref2", 3, b"TT")
    seq2 = _make_seq("seq2", 4, b"TT")
    proc = save_fasta(pairwise=True)
    output = list(proc([(ref1, seq1), (ref2, seq2)], []))
    assert output == [
        ">ref1\n",
        "AA\n",
        ">seq1\n",
        "AA\n",
        ">ref2\n",
        "TT\n",
        ">seq2\n",
        "TT\n",
    ]


def test_save_fasta_preserve_order_skips_ref_for_nonsequential() -> None:
    """Non-sequential IDs omit the reference when preserving order."""
    from postalign.processors.save_fasta import save_fasta

    ref = _make_seq("ref", 1, b"AA")
    seq = _make_seq("seq", 3, b"AA")
    proc = save_fasta(preserve_order=True)
    output = list(proc([(ref, seq)], []))
    assert output == [">seq\n", "AA\n"]
