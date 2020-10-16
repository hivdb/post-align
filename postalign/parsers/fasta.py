from ..models.sequence import Sequence


def load(fp, seqtype):
    header = None
    curseq = ''
    seqid = 0
    for line in fp:
        if line.startswith('>'):
            if header and curseq:
                seqid += 1
                yield Sequence(
                    header=header,
                    seqtext=curseq,
                    seqid=seqid,
                    seqtype=seqtype,
                    skip_invalid=True)
            header = line[1:].strip()
            curseq = ''
        elif line.startswith('#'):
            continue
        else:
            curseq += line.strip()
    if header and curseq:
        seqid += 1
        yield Sequence(
            header=header,
            seqtext=curseq,
            seqid=seqid,
            seqtype=seqtype,
            skip_invalid=True)


def dump(sequences, fp, seqtype):
    pass
