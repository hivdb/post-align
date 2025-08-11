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


def test_paf_load_missing_alignment_yields_empty_pair() -> None:
    """Sequences without matching PAF entries should yield empty texts."""
    from io import StringIO
    from postalign.models import NAPosition, Message
    from postalign.parsers import paf

    paf_stream = StringIO("")
    seqs = StringIO(">1\nAC\n")
    ref = StringIO(">ref\nAC\n")
    messages: list[Message] = []
    refseq, seq = next(
        iter(paf.load(paf_stream, seqs, ref, NAPosition, messages))
    )
    assert refseq.seqtext_as_str == ""
    assert seq.seqtext_as_str == ""


def test_paf_load_inserts_unaligned_region() -> None:
    """Unaligned sequence regions should be inserted with gaps in reference."""
    from io import StringIO
    from postalign.models import NAPosition, Message, PositionFlag
    from postalign.parsers import paf

    paf_text = (
        "1\t6\t0\t2\t+\tref\t4\t0\t2\t2\t2\t60\tcg:Z:2M\n"
        "1\t6\t4\t6\t+\tref\t4\t2\t4\t2\t2\t60\tcg:Z:2M\n"
    )
    paf_stream = StringIO(paf_text)
    seqs = StringIO(">1\nACNNGT\n")
    ref = StringIO(">ref\nACGT\n")
    messages: list[Message] = []
    refseq, seq = next(
        iter(paf.load(paf_stream, seqs, ref, NAPosition, messages))
    )
    assert refseq.seqtext_as_str == "AC--GT"
    assert seq.seqtext_as_str == "ACNNGT"
    assert all(
        pos.flag & PositionFlag.UNALIGNED for pos in seq.seqtext[2:4]
    )


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


def test_paf_load_ref_overlap_skips_alignment() -> None:
    """Fully overlapping reference alignments should be omitted."""
    from io import StringIO
    from types import SimpleNamespace, ModuleType
    from unittest.mock import patch
    import importlib
    import sys
    from postalign.models import NAPosition, Message

    class FakePafRecord(SimpleNamespace):
        @classmethod
        def from_str(cls, line: str) -> "FakePafRecord":
            parts = line.strip().split("\t")
            return cls(
                qname=parts[0],
                qlen=int(parts[1]),
                qstart=int(parts[2]),
                qend=int(parts[3]),
                strand=parts[4],
                tname=parts[5],
                tlen=int(parts[6]),
                tstart=int(parts[7]),
                tend=int(parts[8]),
                mlen=int(parts[9]),
                blen=int(parts[10]),
                mapq=int(parts[11]),
                tags={"cg": SimpleNamespace(value=parts[12].split(":")[-1])},
            )

    fake_pafpy = ModuleType("pafpy")
    fake_pafpy.PafRecord = FakePafRecord  # type: ignore[attr-defined]
    fake_pafpy.Strand = SimpleNamespace(  # type: ignore[attr-defined]
        Reverse="-",
    )

    paf_text = (
        "1\t8\t0\t4\t+\tref\t8\t0\t4\t4\t4\t60\tcg:Z:4M\n"
        "1\t8\t4\t8\t+\tref\t8\t0\t4\t4\t4\t60\tcg:Z:4M\n"
    )
    paf_stream = StringIO(paf_text)
    seqs = StringIO(">1\nAAAAAAAA\n")
    ref = StringIO(">ref\nAAAAAAAA\n")
    messages: list[Message] = []
    with patch.dict(sys.modules, {"pafpy": fake_pafpy}):
        import postalign.parsers.paf as paf_module
        paf_module = importlib.reload(paf_module)
        refseq, _ = next(
            iter(paf_module.load(paf_stream, seqs, ref, NAPosition, messages))
        )
        assert any(
            "REF has already been aligned" in m.message for m in messages
        )
        assert ";" not in refseq.modifiers.last_modifier.text
        paf_module = importlib.reload(paf_module)


def test_paf_load_masks_reused_unaligned_positions() -> None:
    """Unaligned segments reusing aligned positions become gaps."""
    from io import StringIO
    from types import SimpleNamespace, ModuleType
    from unittest.mock import patch
    import importlib
    import sys
    from postalign.models import NAPosition, Sequence, Message

    ref_text = NAPosition.init_from_bytes(b"AAAAAAAA")
    refseq = Sequence(
        header="ref",
        description="",
        seqtext=ref_text,
        seqid=0,
        seqtype=NAPosition,
        abs_seqstart=0,
    )
    seq_text = NAPosition.init_from_bytes(b"AAAATTTT")
    for idx, pos in enumerate(seq_text[4:], start=1):
        pos.pos = idx
    seq = Sequence(
        header="seq",
        description="",
        seqtext=seq_text,
        seqid=1,
        seqtype=NAPosition,
        abs_seqstart=0,
    )

    fake_pafpy = ModuleType("pafpy")
    fake_pafpy.PafRecord = SimpleNamespace(  # type: ignore[attr-defined]
        from_str=lambda line: SimpleNamespace(
            qname="1",
            qlen=8,
            qstart=0,
            qend=4,
            strand="+",
            tname="ref",
            tlen=8,
            tstart=0,
            tend=4,
            mlen=4,
            blen=4,
            mapq=60,
            tags={"cg": SimpleNamespace(value="4M")},
        )
    )
    fake_pafpy.Strand = SimpleNamespace(  # type: ignore[attr-defined]
        Reverse="-",
    )

    paf_stream = StringIO(
        "1\t8\t0\t4\t+\tref\t8\t0\t4\t4\t4\t60\tcg:Z:4M\n"
    )
    messages: list[Message] = []
    with patch.dict(sys.modules, {"pafpy": fake_pafpy}), patch(
        "postalign.parsers.paf.fasta.load",
        side_effect=[iter([refseq]), iter([seq])],
    ):
        import postalign.parsers.paf as paf_module
        paf_module = importlib.reload(paf_module)
        _, seq_out = next(
            iter(
                paf_module.load(
                    paf_stream,
                    StringIO(""),
                    StringIO(""),
                    NAPosition,
                    messages,
                )
            )
        )
        assert seq_out.seqtext_as_str.endswith("----")
        paf_module = importlib.reload(paf_module)
