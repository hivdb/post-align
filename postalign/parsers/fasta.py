from typing import Type, TextIO, Iterable, Generator

from ..models import Sequence, Position

GAP_CHARS = b'.-'


def load(
    fp: TextIO,
    seqtype: Type[Position],
    *,
    remove_gaps: bool = False
) -> Generator[Sequence, None, None]:
    header: str = ''
    curseq: bytearray = bytearray()
    seqid: int = 0

    def make_seq() -> Sequence:
        nonlocal seqid
        seqid += 1
        if remove_gaps:
            for gap in GAP_CHARS:
                curseq.replace(bytes([gap]), b'')

        headerdesc = header.split(' ', 1)
        description = ''
        if len(headerdesc) == 2:
            description = headerdesc[1]
        return Sequence(
            header=headerdesc[0],
            description=description,
            seqtext=seqtype.init_from_bytes(curseq),
            seqid=seqid,
            seqtype=seqtype,
            abs_seqstart=0,
            skip_invalid=True)

    for line in fp:
        if line.startswith('>'):
            if header:
                yield make_seq()
            header = line[1:].strip()
            curseq = bytearray()
        elif line.startswith('#'):
            continue
        else:
            curseq.extend(bytes(line.strip(), 'ASCII', 'ignore'))
    if header:
        yield make_seq()


def dump(
    sequences: Iterable[Sequence],
    fp: TextIO,
    seqtype: str
) -> None:
    pass
