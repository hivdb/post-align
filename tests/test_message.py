"""Tests for the :mod:`postalign.models.message` module."""

from __future__ import annotations


def test_message_string_and_dict() -> None:
    """Ensure message string representation and dictionary conversion work."""
    from postalign.models.message import Message, MessageLevel

    message = Message(1, MessageLevel.WARNING, "issue")
    assert str(message) == "[WARNING] 1:issue"
    assert repr(message) == "<Message [WARNING] 1:issue>"
    assert message.to_dict() == {"level": "WARNING", "message": "issue"}
