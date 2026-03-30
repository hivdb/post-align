from postalign_rs import NAPosition, NAPositionList, enumerate_seq_pos

GAP_CHAR: int = ord(b'-')
GAP_CHARS: tuple[int, ...] = tuple(b'-.')

__all__ = [
    'GAP_CHAR',
    'GAP_CHARS',
    'NAPosition',
    'NAPositionList',
    'enumerate_seq_pos',
]
