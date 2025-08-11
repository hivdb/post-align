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


def test_insert_unaligned_region_zero_sizes_noop() -> None:
    """Zero-sized unaligned regions should be ignored."""
    from postalign.parsers.paf import insert_unaligned_region
    from postalign.models import NAPosition

    reftext = NAPosition.init_from_bytes(b"AB")
    seqtext = NAPosition.init_from_bytes(b"AB")
    orig = seqtext[:]
    insert_unaligned_region(
        reftext,
        seqtext,
        orig,
        NAPosition,
        align1_ref_end=1,
        align1_seq_end=1,
        align2_ref_start=1,
        align2_seq_start=1,
    )
    assert NAPosition.as_str(reftext) == "AB"
    assert NAPosition.as_str(seqtext) == "AB"


def test_insert_unaligned_region_defaults_to_alignment2() -> None:
    """Default placement should favor the second alignment."""
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
    )
    assert NAPosition.as_str(reftext) == "ABCD--"
    assert NAPosition.as_str(seqtext) == "ABCDXX"


def test_insert_unaligned_region_seq_only() -> None:
    """Reference gaps should be inserted when only sequence extends."""
    from postalign.parsers.paf import insert_unaligned_region
    from postalign.models import NAPosition

    reftext = NAPosition.init_from_bytes(b"AB")
    seqtext = NAPosition.init_from_bytes(b"AB")
    orig = NAPosition.init_from_bytes(b"AXB")
    insert_unaligned_region(
        reftext,
        seqtext,
        orig,
        NAPosition,
        align1_ref_end=1,
        align1_seq_end=1,
        align2_ref_start=1,
        align2_seq_start=2,
    )
    assert NAPosition.as_str(reftext) == "A-B"
    assert NAPosition.as_str(seqtext) == "AXB"


def test_insert_unaligned_region_ref_only_noop() -> None:
    """Reference-only gaps should leave texts unchanged."""
    from postalign.parsers.paf import insert_unaligned_region
    from postalign.models import NAPosition

    reftext = NAPosition.init_from_bytes(b"AB")
    seqtext = NAPosition.init_from_bytes(b"AB")
    orig = seqtext[:]
    insert_unaligned_region(
        reftext,
        seqtext,
        orig,
        NAPosition,
        align1_ref_end=1,
        align1_seq_end=1,
        align2_ref_start=3,
        align2_seq_start=1,
    )
    assert NAPosition.as_str(reftext) == "AB"
    assert NAPosition.as_str(seqtext) == "AB"


def test_paf_load_no_alignment_returns_error() -> None:
    """Sequences without matching PAF entries should yield empty texts."""
    from io import StringIO
    from unittest.mock import patch
    import importlib
    import sys
    import pafpy  # type: ignore[import-untyped]
    from postalign.models import NAPosition, Message

    paf_stream = StringIO("")
    seqs = StringIO(">1\nAC\n")
    ref = StringIO(">ref\nAC\n")
    messages: list[Message] = []
    with patch.dict(sys.modules, {"pafpy": pafpy}):
        import postalign.parsers.paf as paf_module
        paf_module = importlib.reload(paf_module)
        refseq, seq = list(
            paf_module.load(paf_stream, seqs, ref, NAPosition, messages)
        )[0]
        assert refseq.seqtext_as_str == ""
        assert seq.seqtext_as_str == ""
        paf_module = importlib.reload(paf_module)


def test_paf_load_skips_reverse_strand() -> None:
    """Reverse-strand alignments should be ignored during loading."""
    from io import StringIO
    from types import SimpleNamespace, ModuleType
    from unittest.mock import patch
    import importlib
    import sys
    from postalign.models import NAPosition, Message

    fake_pafpy = ModuleType("pafpy")

    class FakePafRecord:
        def __init__(
            self,
            qname: str,
            tstart: int,
            tend: int,
            qstart: int,
            qend: int,
            strand: str,
            cigar: str,
        ) -> None:
            self.qname = qname
            self.tstart = tstart
            self.tend = tend
            self.qstart = qstart
            self.qend = qend
            self.strand = strand
            self.tags = {"cg": SimpleNamespace(value=cigar)}

        @classmethod
        def from_str(cls, line: str) -> "FakePafRecord":
            parts = line.strip().split("\t")
            return cls(
                parts[0],
                int(parts[7]),
                int(parts[8]),
                int(parts[2]),
                int(parts[3]),
                parts[4],
                parts[12].split(":")[-1],
            )

    fake_pafpy.PafRecord = FakePafRecord  # type: ignore[attr-defined]
    fake_pafpy.Strand = SimpleNamespace(  # type: ignore[attr-defined]
        Reverse="-"
    )

    paf_stream = StringIO(
        "1\t4\t0\t4\t-\tref\t4\t0\t4\t4\t4\t60\tcg:Z:4M\n"
    )
    seqs = StringIO(">1\nACGT\n")
    ref = StringIO(">ref\nACGT\n")
    messages: list[Message] = []
    with patch.dict(sys.modules, {"pafpy": fake_pafpy}):
        import postalign.parsers.paf as paf_module
        paf_module = importlib.reload(paf_module)
        refseq, seq = next(
            iter(paf_module.load(paf_stream, seqs, ref, NAPosition, messages))
        )
        assert seq.seqtext_as_str == ""
        paf_module = importlib.reload(paf_module)


def test_paf_load_overlapping_seq_range_warns() -> None:
    """Overlapping sequence ranges should emit a warning."""
    from io import StringIO
    from types import SimpleNamespace, ModuleType
    from unittest.mock import patch
    import importlib
    import sys
    from postalign.models import NAPosition, Message, MessageLevel

    fake_pafpy = ModuleType("pafpy")

    class FakePafRecord:
        def __init__(
            self,
            qname: str,
            tstart: int,
            tend: int,
            qstart: int,
            qend: int,
            strand: str,
            cigar: str,
        ) -> None:
            self.qname = qname
            self.tstart = tstart
            self.tend = tend
            self.qstart = qstart
            self.qend = qend
            self.strand = strand
            self.tags = {"cg": SimpleNamespace(value=cigar)}

        @classmethod
        def from_str(cls, line: str) -> "FakePafRecord":
            parts = line.strip().split("\t")
            return cls(
                parts[0],
                int(parts[7]),
                int(parts[8]),
                int(parts[2]),
                int(parts[3]),
                parts[4],
                parts[12].split(":")[-1],
            )

    fake_pafpy.PafRecord = FakePafRecord  # type: ignore[attr-defined]
    fake_pafpy.Strand = SimpleNamespace(  # type: ignore[attr-defined]
        Reverse="-"
    )

    paf_text = (
        "1\t8\t0\t4\t+\tref\t8\t0\t4\t4\t4\t60\tcg:Z:4M\n"
        "1\t8\t2\t6\t+\tref\t8\t4\t8\t4\t4\t60\tcg:Z:4M\n"
    )
    paf_stream = StringIO(paf_text)
    seqs = StringIO(">1\nAAAAAAAA\n")
    ref = StringIO(">ref\nAAAAAAAA\n")
    messages: list[Message] = []
    with patch.dict(sys.modules, {"pafpy": fake_pafpy}):
        import postalign.parsers.paf as paf_module
        paf_module = importlib.reload(paf_module)
        list(paf_module.load(paf_stream, seqs, ref, NAPosition, messages))
        assert any(
            "SEQ has already been aligned" in m.message for m in messages
        )
        assert all(m.level == MessageLevel.WARNING for m in messages)
        paf_module = importlib.reload(paf_module)


def test_paf_load_ref_overlap_shrinks_cigar() -> None:
    """Reference overlaps should shrink the CIGAR and warn the user."""
    from io import StringIO
    from types import SimpleNamespace, ModuleType
    from unittest.mock import patch
    import importlib
    import sys
    from postalign.models import NAPosition, Message

    fake_pafpy = ModuleType("pafpy")

    class FakePafRecord:
        def __init__(
            self,
            qname: str,
            tstart: int,
            tend: int,
            qstart: int,
            qend: int,
            strand: str,
            cigar: str,
        ) -> None:
            self.qname = qname
            self.tstart = tstart
            self.tend = tend
            self.qstart = qstart
            self.qend = qend
            self.strand = strand
            self.tags = {"cg": SimpleNamespace(value=cigar)}

        @classmethod
        def from_str(cls, line: str) -> "FakePafRecord":
            parts = line.strip().split("\t")
            return cls(
                parts[0],
                int(parts[7]),
                int(parts[8]),
                int(parts[2]),
                int(parts[3]),
                parts[4],
                parts[12].split(":")[-1],
            )

    fake_pafpy.PafRecord = FakePafRecord  # type: ignore[attr-defined]
    fake_pafpy.Strand = SimpleNamespace(  # type: ignore[attr-defined]
        Reverse="-"
    )

    paf_text = (
        "1\t8\t0\t4\t+\tref\t8\t0\t4\t4\t4\t60\tcg:Z:4M\n"
        "1\t8\t4\t8\t+\tref\t8\t2\t6\t4\t4\t60\tcg:Z:4M\n"
    )
    paf_stream = StringIO(paf_text)
    seqs = StringIO(">1\nAAAAAAAA\n")
    ref = StringIO(">ref\nAAAAAAAA\n")
    messages: list[Message] = []
    with patch.dict(sys.modules, {"pafpy": fake_pafpy}):
        import postalign.parsers.paf as paf_module
        paf_module = importlib.reload(paf_module)
        refseq, seq = next(
            iter(paf_module.load(paf_stream, seqs, ref, NAPosition, messages))
        )
        assert any(
            "REF has already been aligned" in m.message for m in messages
        )
        # The first alignment should be shrunk to length 2
        assert "0,2,2M" in refseq.modifiers.last_modifier.text
        paf_module = importlib.reload(paf_module)
