import re
from itertools import groupby
from ..models.sequence import PositionalSeqStr

CIGAR_PATTERN = re.compile(r'(\d+)([MNDI])')


def get_cigar(ref, query):
    cigar = []
    for r, q in zip(ref, query):
        if r == q == '-':
            continue
        if r == '-':
            cigar.append(('I', 1))
        elif q == '-':
            cigar.append(('D', 1))
        else:
            cigar.append(('M', 1))
    compressed = []
    for type, counts in groupby(cigar, lambda c: c[0]):
        total = sum(c for _, c in counts)
        compressed.append('{}{}'.format(total, type))
    return ''.join(compressed)


class CIGAR:

    def __init__(self, ref_start, seq_start, cigar_string):
        self.ref_start = ref_start
        self.seq_start = seq_start
        self.cigar_string = cigar_string
        self.cigar_tuple = [
            (int(num), op)
            for num, op in CIGAR_PATTERN.findall(cigar_string)
        ]

    def get_alignment(self, refseq, seq):
        aligned_refseq = refseq[self.ref_start:]
        aligned_seq = seq[self.seq_start:]
        offset = 0
        for num, op in self.cigar_tuple:
            if op == 'M':
                offset += num
            elif op in ('D', 'N'):
                aligned_seq = (aligned_seq[:offset] +
                               PositionalSeqStr('-' * num) +
                               aligned_seq[offset:])
                offset += num
            elif op == 'I':
                aligned_refseq = (aligned_refseq[:offset] +
                                  PositionalSeqStr('-' * num) +
                                  aligned_refseq[offset:])
                offset += num
        aligned_seq = (PositionalSeqStr('-' * (self.ref_start)) +
                       aligned_seq[:offset] +
                       PositionalSeqStr('-' * (len(aligned_refseq) - offset)))
        aligned_refseq = refseq[:self.ref_start] + aligned_refseq
        if len(aligned_refseq) != len(aligned_seq):
            raise ValueError(
                'Unmatched alignment length: {!r} and {!r}'
                .format(aligned_refseq, aligned_seq)
            )
        return aligned_refseq, aligned_seq

    def __repr__(self):
        return ('<CIGAR {!r} ref_start={!r} seq_start={!r}>'
                .format(self.cigar_string,
                        self.ref_start,
                        self.seq_start))
