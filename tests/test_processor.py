"""Tests for :mod:`postalign.processor`."""

from __future__ import annotations

from unittest.mock import MagicMock


def test_output_processor_decorator() -> None:
    """Output processors should delegate calls and set flags."""
    from postalign.processor import output_processor
    from postalign.models import Message

    func = MagicMock(return_value=['out'])
    proc = output_processor('cmd')(func)
    messages: list[Message] = []
    result = list(proc([], messages))
    func.assert_called_once_with([], messages)
    assert proc.is_output_command
    assert result == ['out']


def test_intermediate_processor_decorator() -> None:
    """Intermediate processors propagate results and disable output flag."""
    from postalign.processor import intermediate_processor
    from postalign.models import Message

    func = MagicMock(return_value=['mid'])
    proc = intermediate_processor('cmd')(func)
    messages: list[Message] = []
    result = list(proc([], messages))
    func.assert_called_once_with([], messages)
    assert not proc.is_output_command
    assert result == ['mid']
