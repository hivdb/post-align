"""Tests for :mod:`postalign.cli` helper functions."""

from __future__ import annotations

import pytest
import typer
from typing import Callable
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, patch


def _noop_result_callback(
    *_: object,
    **__: object,
) -> Callable[[Callable[..., object]], Callable[..., object]]:
    """Return a decorator that leaves the function unchanged."""

    def decorator(func: Callable[..., object]) -> Callable[..., object]:
        return func

    return decorator


def test_call_processors_runs_pipeline() -> None:
    """Processors should run sequentially and yield final output."""
    with patch.object(
        typer.Typer,
        "result_callback",
        _noop_result_callback,
        create=True,
    ):
        from postalign.cli import call_processors
        from postalign.processor import Processor
        from postalign.models import Message
        from postalign.models.sequence import Sequence, RefSeqPair

    seq_a = MagicMock(spec=Sequence)
    seq_b = MagicMock(spec=Sequence)
    iterator: list[RefSeqPair] = [(seq_a, seq_b)]
    mid_seq_a = MagicMock(spec=Sequence)
    mid_seq_b = MagicMock(spec=Sequence)
    mid_iterator: list[RefSeqPair] = [(mid_seq_a, mid_seq_b)]
    first_func = MagicMock(return_value=mid_iterator)
    last_func = MagicMock(return_value=['out'])
    first = Processor('first', False, first_func)
    last = Processor('last', True, last_func)
    messages: list[Message] = []
    result = list(call_processors([first, last], iterator, messages))
    first_func.assert_called_once_with(iterator, messages)
    last_func.assert_called_once_with(mid_iterator, messages)
    assert result == ['out']


def test_check_processors_errors() -> None:
    """`check_processors` should validate pipeline configuration."""
    with patch.object(
        typer.Typer,
        "result_callback",
        _noop_result_callback,
        create=True,
    ):
        from postalign.cli import check_processors
        from postalign.processor import Processor

    with pytest.raises(typer.BadParameter):
        check_processors([])
    mid = Processor('mid', False, MagicMock())
    with pytest.raises(typer.BadParameter):
        check_processors([mid])
    p1 = Processor('p1', True, MagicMock())
    p2 = Processor('p2', True, MagicMock())
    with pytest.raises(typer.BadParameter):
        check_processors([p1, p2])


def test_reference_callback_missing_name() -> None:
    """Missing parameter metadata should raise :class:`BadParameter`."""
    with patch.object(
        typer.Typer,
        "result_callback",
        _noop_result_callback,
        create=True,
    ):
        from postalign.cli import AlignmentFormat, reference_callback

    ctx = MagicMock()
    ctx.params = {'alignment_format': AlignmentFormat.MSA}
    param = MagicMock()
    param.name = None
    with pytest.raises(typer.BadParameter):
        reference_callback(ctx, param, 'ref')


def test_seqs_prior_alignment_callback_requires_file() -> None:
    """PAF format requires a provided file handle."""
    with patch.object(
        typer.Typer,
        "result_callback",
        _noop_result_callback,
        create=True,
    ):
        from postalign.cli import (
            AlignmentFormat,
            seqs_prior_alignment_callback,
        )

    ctx = MagicMock()
    ctx.params = {
        'alignment_format': AlignmentFormat.PAF,
    }
    param = MagicMock()
    param.name = 'seqs_prior_alignment'
    with pytest.raises(typer.BadParameter):
        seqs_prior_alignment_callback(ctx, param, None)


def test_reference_callback_reads_msa_header(tmp_path: Path) -> None:
    """MSA references should return the parsed header text."""
    with patch.object(
        typer.Typer,
        "result_callback",
        _noop_result_callback,
        create=True,
    ):
        from postalign.cli import AlignmentFormat, reference_callback

    ctx = MagicMock()
    ctx.params = {"alignment_format": AlignmentFormat.MSA}
    param = MagicMock()
    param.name = "reference"
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">ref\nACG\n")
    result = reference_callback(ctx, param, str(fasta))
    assert result == "ref"


def test_seqs_prior_alignment_callback_ignores_for_non_paf() -> None:
    """Non-PAF formats should ignore provided files with a warning."""
    with patch.object(
        typer.Typer,
        "result_callback",
        _noop_result_callback,
        create=True,
    ):
        from postalign.cli import (
            AlignmentFormat,
            seqs_prior_alignment_callback,
        )

    ctx = MagicMock()
    ctx.params = {"alignment_format": AlignmentFormat.MSA}
    param = MagicMock()
    param.name = "seqs_prior_alignment"
    buf = StringIO()
    with patch("sys.stderr", new=StringIO()) as err:
        result = seqs_prior_alignment_callback(ctx, param, buf)
    assert result is None
    assert "ignore -p/--seqs-prior-alignment" in err.getvalue()


def test_reference_callback_requires_file_for_non_msa() -> None:
    """Non-MSA formats should enforce file-based references."""
    with patch.object(
        typer.Typer,
        "result_callback",
        _noop_result_callback,
        create=True,
    ):
        from postalign.cli import AlignmentFormat, reference_callback
    ctx = MagicMock()
    ctx.params = {"alignment_format": AlignmentFormat.PAF}
    param = MagicMock()
    param.name = "reference"
    with patch("builtins.open", side_effect=OSError), pytest.raises(
        typer.BadParameter
    ):
        reference_callback(ctx, param, "ref.fa")
