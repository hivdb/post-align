"""Tests for the :mod:`postalign.models.message` module."""

from postalign.models.message import Message, MessageLevel


def test_message_string_and_dict() -> None:
    """Ensure message string representation and dictionary conversion work."""
    message = Message(1, MessageLevel.WARNING, "issue")
    assert str(message) == "[WARNING] 1:issue"
    assert repr(message) == "<Message [WARNING] 1:issue>"
    assert message.to_dict() == {"level": "WARNING", "message": "issue"}
