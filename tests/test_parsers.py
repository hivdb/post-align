"""Tests for parser utilities."""

from __future__ import annotations

from io import StringIO
from unittest.mock import MagicMock, patch

import pytest


def test_fasta_load_removes_gaps() -> None:
    """FASTA loader should strip gaps and parse descriptions."""
    from postalign.parsers import fasta
    from postalign.models import NAPosition

    fp = StringIO(">s1 first\nAC-G\n#comment\n>s2\nTT--AA\n")
    seqs = list(fasta.load(fp, NAPosition, remove_gaps=True))
    assert seqs[0].header == "s1"
    assert seqs[0].description == "first"
    assert seqs[0].seqtext_as_str == "ACG"
    assert seqs[1].seqid == 2
    assert seqs[1].seqtext_as_str == "TTAA"


def test_msa_load_yields_pairs() -> None:
    """MSA loader should yield reference/query pairs."""
    from postalign.parsers import msa
    from postalign.models import NAPosition

    fp = StringIO(">ref\nACG\n>query\nAC-\n")
    pairs = list(msa.load(fp, reference="ref", seqtype=NAPosition))
    ref, seq = pairs[0]
    assert ref.header == "ref"
    assert seq.header == "query"


def test_msa_load_missing_reference() -> None:
    """Missing reference sequence should raise :class:`BadParameter`."""
    from postalign.parsers import msa
    from postalign.models import NAPosition
    import typer

    fp = StringIO(">ref\nACG\n>query\nAC-\n")
    with pytest.raises(typer.BadParameter):
        list(msa.load(fp, reference="missing", seqtype=NAPosition))


def test_minimap2_load_errors_on_failure() -> None:
    """Non-zero minimap2 exit should surface as :class:`BadParameter`."""
    from postalign.parsers import minimap2
    from postalign.models import NAPosition, Message
    import typer

    fastafp = StringIO(">1\nACG\n")
    ref = StringIO(">ref\nACG\n")
    messages: list[Message] = []

    proc = MagicMock()
    proc.communicate.return_value = ("", "boom")
    proc.returncode = 1
    with patch("postalign.parsers.minimap2.Popen", return_value=proc):
        with pytest.raises(typer.BadParameter):
            list(minimap2.load(fastafp, ref, NAPosition, messages))


def test_minimap2_load_handles_timeout() -> None:
    """Timeouts should kill the process and still return records."""
    from postalign.parsers import minimap2
    from postalign.models import NAPosition, Message
    from subprocess import TimeoutExpired

    fastafp = StringIO(">1\nACG\n")
    ref = StringIO(">ref\nACG\n")
    messages: list[Message] = []

    proc = MagicMock()
    proc.communicate.side_effect = [
        TimeoutExpired(cmd="mm", timeout=1),
        ("", ""),
    ]
    proc.kill = MagicMock()
    proc.returncode = 0
    with patch("postalign.parsers.minimap2.Popen", return_value=proc):
        list(minimap2.load(fastafp, ref, NAPosition, messages))
    proc.kill.assert_called_once()


def test_msa_load_uses_first_sequence_as_reference() -> None:
    """MSA loader defaults to first sequence when reference missing."""
    from postalign.parsers import msa
    from postalign.models import NAPosition

    fp = StringIO(">ref\nACG\n>query\nAC-\n")
    pairs = list(msa.load(fp, reference="", seqtype=NAPosition))
    refseq, seq = pairs[0]
    assert refseq.header == "ref"
    assert seq.header == "query"
