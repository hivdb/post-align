"""Tests for the :class:`postalign.models.sequence.Sequence` class."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:  # pragma: no cover - used only for type checking
    from postalign.models import Sequence


def _make_sequence(text: bytes, seqid: int = 1) -> Sequence:
    """Utility to build a :class:`Sequence` from raw bytes."""
    from postalign.models import NAPosition, Sequence
    from postalign.models._sequence import SKIP_VALIDATION

    nas = NAPosition.init_from_bytes(text)
    return Sequence(
        header="h",
        description="d",
        seqtext=nas,
        seqid=seqid,
        seqtype=NAPosition,
        abs_seqstart=1,
        skip_invalid=SKIP_VALIDATION,
    )


def test_sequence_headerdesc_and_modifiers() -> None:
    """Header properties should include description and modifiers."""
    seq = _make_sequence(b"AT")
    assert seq.headerdesc == "h d"
    new_seq = seq.push_seqtext(seq.seqtext + seq.seqtext, "dup", len(seq))
    assert new_seq.abs_seqstart == 3
    assert new_seq.header_with_modifiers == "h MOD::1:dup"


def test_sequence_add_mismatched_seqid_raises() -> None:
    """Adding sequences with different IDs should raise ``ValueError``."""
    seq1 = _make_sequence(b"AT", seqid=1)
    seq2 = _make_sequence(b"AT", seqid=2)
    with pytest.raises(ValueError):
        _ = seq1 + seq2


def test_sequence_getitem_step_slice_error() -> None:
    """Slicing with a step is unsupported and should raise ``ValueError``."""
    seq = _make_sequence(b"AT")
    with pytest.raises(ValueError):
        _ = seq[::2]
