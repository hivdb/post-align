from typing import TextIO, Iterable, Generator

from ..models.sequence import Sequence, PositionalSeqStr

GAP_CHARS = '.-'


def load(
    fp: TextIO,
    seqtype: str,
    *,
    remove_gaps: bool = False
) -> Generator[Sequence, None, None]:
    header: str = ''
    curseq: str = ''
    seqid: int = 0

    def make_seq() -> Sequence:
        nonlocal seqid
        seqid += 1
        if remove_gaps:
            for gap in GAP_CHARS:
                curseq.replace(gap, '')

        headerdesc = header.split(' ', 1)
        description = ''
        if len(headerdesc) == 2:
            description = headerdesc[1]
        return Sequence(
            header=headerdesc[0],
            description=description,
            seqtext=PositionalSeqStr.init_from_nastring(curseq),
            seqid=seqid,
            seqtype=seqtype,
            abs_seqstart=0,
            skip_invalid=True)

    for line in fp:
        if line.startswith('>'):
            if header and curseq:
                yield make_seq()
            header = line[1:].strip()
            curseq = ''
        elif line.startswith('#'):
            continue
        else:
            curseq += line.strip()
    if header and curseq:
        yield make_seq()


def dump(
    sequences: Iterable[Sequence],
    fp: TextIO,
    seqtype: str
) -> None:
    pass
