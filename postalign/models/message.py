from enum import Enum
from typing import Dict


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
        return '[{}] {}:{}'.format(self.level.name, self.seqid, self.message)

    def __repr__(self) -> str:
        return '<Message {}>'.format(self)

    def to_dict(self) -> Dict[str, str]:
        return {
            'level': self.level.name,
            'message': self.message
        }
