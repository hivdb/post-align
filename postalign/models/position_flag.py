from enum import IntFlag


class PositionFlag(IntFlag):
    UNALIGNED = 0x10
    TRIM_BY_SEQ = 0x01
    NONE = 0x00
