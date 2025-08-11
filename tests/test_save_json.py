"""Tests for JSON saving utilities."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest


def test_gene_range_tuples_callback_parses_ranges() -> None:
    """Callback should group gene ranges into tuples."""
    from postalign.processors.save_json import gene_range_tuples_callback

    ctx = MagicMock()
    param = MagicMock()
    result = gene_range_tuples_callback(ctx, param, ("G", "1", "4"))
    assert result == [("G", [(1, 4)])]


def test_gene_range_tuples_callback_missing_value() -> None:
    """Odd range values should raise :class:`BadParameter`."""
    from postalign.processors.save_json import gene_range_tuples_callback
    import typer

    ctx = MagicMock()
    param = MagicMock()
    with pytest.raises(typer.BadParameter):
        gene_range_tuples_callback(ctx, param, ("G", "1"))


def test_gene_range_tuples_callback_refstart_too_small() -> None:
    """Reference start below 1 should raise :class:`BadParameter`."""
    from postalign.processors.save_json import gene_range_tuples_callback
    import typer

    ctx = MagicMock()
    param = MagicMock()
    with pytest.raises(typer.BadParameter):
        gene_range_tuples_callback(ctx, param, ("G", "0", "4"))


def test_gene_range_tuples_callback_refend_too_small() -> None:
    """Reference end must allow at least one codon."""
    from postalign.processors.save_json import gene_range_tuples_callback
    import typer

    ctx = MagicMock()
    param = MagicMock()
    with pytest.raises(typer.BadParameter):
        gene_range_tuples_callback(ctx, param, ("G", "2", "3"))


def test_gene_range_tuples_callback_multiple_ranges() -> None:
    """Callback should support multiple ranges for a gene."""
    from postalign.processors.save_json import gene_range_tuples_callback

    ctx = MagicMock()
    param = MagicMock()
    result = gene_range_tuples_callback(
        ctx, param, ("G", "1", "4", "7", "10")
    )
    assert result == [("G", [(1, 4), (7, 10)])]


def test_gene_range_tuples_callback_two_genes() -> None:
    """Multiple genes should yield separate tuples."""
    from postalign.processors.save_json import gene_range_tuples_callback

    ctx = MagicMock()
    param = MagicMock()
    result = gene_range_tuples_callback(
        ctx, param, ("G1", "1", "4", "G2", "5", "8")
    )
    assert result == [("G1", [(1, 4)]), ("G2", [(5, 8)])]


def test_gene_range_tuples_callback_empty_input() -> None:
    """Callback should handle an empty argument list."""
    from postalign.processors.save_json import gene_range_tuples_callback

    ctx = MagicMock()
    param = MagicMock()
    assert gene_range_tuples_callback(ctx, param, ()) == []


def test_save_json_includes_frameshift() -> None:
    """Processor output should list frameshift events."""
    from unittest.mock import patch
    from postalign.processors.save_json import save_json
    from postalign.models import Sequence, NAPosition, Message

    refcd = NAPosition.init_from_bytes(b"AAA")
    seqcd = refcd + [NAPosition.init_from_bytes(b"T")[0]]
    gene_codons = [("G", [refcd], [seqcd])]
    with patch(
        "postalign.processors.save_json.group_by_gene_codons",
        return_value=gene_codons,
    ), patch(
        "postalign.processors.save_json.find_codon_trim_slice",
        return_value=slice(0, 1),
    ), patch(
        "postalign.processors.save_json.translate_codon",
        return_value=b"M",
    ):
        processor = save_json([("G", [(1, 4)])])
        refseq = Sequence(
            header=">ref",
            description="",
            seqtext=refcd,
            seqid=1,
            seqtype=NAPosition,
            abs_seqstart=0,
        )
        seq = Sequence(
            header=">seq",
            description="",
            seqtext=seqcd,
            seqid=1,
            seqtype=NAPosition,
            abs_seqstart=0,
        )
        messages: list[Message] = []
        output = "".join(processor([(refseq, seq)], messages))
        assert '"FrameShifts":' in output


def test_save_json_reports_unaligned_gene() -> None:
    """Genes without codons should record an error message."""
    from unittest.mock import patch
    from postalign.processors.save_json import save_json
    from postalign.models import Sequence, NAPosition, Message

    with patch(
        "postalign.processors.save_json.group_by_gene_codons",
        return_value=[("G", [], [])],
    ), patch(
        "postalign.processors.save_json.find_codon_trim_slice",
        return_value=slice(0, 0),
    ):
        processor = save_json([("G", [(1, 4)])])
        refseq = Sequence(
            header=">ref",
            description="",
            seqtext=[],
            seqid=1,
            seqtype=NAPosition,
            abs_seqstart=0,
        )
        seq = Sequence(
            header=">seq",
            description="",
            seqtext=[],
            seqid=1,
            seqtype=NAPosition,
            abs_seqstart=0,
        )
        messages: list[Message] = []
        output = "".join(processor([(refseq, seq)], messages))
        assert '"Error": "Sequence is not aligned"' in output
