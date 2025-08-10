from enum import Enum


class MessageLevel(Enum):
    INFO = 0
    WARNING = 1
    ERROR = 2


class Message:

    def __init__(self, seqid: int, level: MessageLevel, message: str):
        self.seqid = seqid
        self.level = level
        self.message = message

    def __str__(self) -> str:
        return f'[{self.level.name}] {self.seqid}:{self.message}'

    def __repr__(self) -> str:
        return f'<Message {self}>'

    def to_dict(self) -> dict[str, str]:
        return {
            'level': self.level.name,
            'message': self.message
        }
