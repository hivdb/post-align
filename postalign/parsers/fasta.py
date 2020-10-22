from ..models.sequence import Sequence

GAP_CHARS = '.-'


def load(fp, seqtype, *, remove_gaps=False):
    header = None
    curseq = ''
    seqid = 0

    def make_seq():
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
            seqtext=curseq,
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


def dump(sequences, fp, seqtype):
    pass
