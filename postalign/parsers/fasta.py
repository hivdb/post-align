from ..models.sequence import Sequence

GAP_CHARS = '.-'


def load(fp, seqtype, *, remove_gaps=False):
    header = None
    curseq = ''
    seqid = 0
    for line in fp:
        if line.startswith('>'):
            if header and curseq:
                seqid += 1
                if remove_gaps:
                    for gap in GAP_CHARS:
                        curseq.replace(gap, '')

                yield Sequence(
                    header=header,
                    seqtext=curseq,
                    seqid=seqid,
                    seqtype=seqtype,
                    abs_seqstart=0,
                    skip_invalid=True)
            header = line[1:].strip()
            curseq = ''
        elif line.startswith('#'):
            continue
        else:
            curseq += line.strip()
    if header and curseq:
        seqid += 1
        if remove_gaps:
            for gap in GAP_CHARS:
                curseq.replace(gap, '')

        yield Sequence(
            header=header,
            seqtext=curseq,
            seqid=seqid,
            seqtype=seqtype,
            abs_seqstart=0,
            skip_invalid=True)


def dump(sequences, fp, seqtype):
    pass
