"""Tests for the :mod:`postalign.models.message` module."""

from __future__ import annotations


def test_message_string_and_dict() -> None:
    """Ensure message string representation and dictionary conversion work."""
    from postalign.models.message import Message, MessageLevel

    message = Message(1, MessageLevel.WARNING, "issue")
    assert str(message) == "[WARNING] 1:issue"
    assert repr(message) == "<Message [WARNING] 1:issue>"
    assert message.to_dict() == {"level": "WARNING", "message": "issue"}


def test_message_level_enum_values() -> None:
    """Enum members should expose their integer levels."""
    from postalign.models.message import MessageLevel

    assert MessageLevel.INFO.value == 0
    assert MessageLevel.ERROR.value == 2
