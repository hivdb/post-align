import re

from ..cli import cli

LEFT_TRIM_PATTERN = re.compile(r'^[.-]+')
RIGHT_TRIM_PATTERN = re.compile(r'[.-]+$')


def find_trim_slice(seq):
    left_trim = None
    match = LEFT_TRIM_PATTERN.search(seq.seqtext)
    if match:
        left_trim = match.span()[1]
    right_trim = None
    match = RIGHT_TRIM_PATTERN.search(seq.seqtext)
    if match:
        right_trim = match.span()[0]
    return slice(left_trim, right_trim)


@cli.command('trim-by-ref')
def trim_by_ref():
    """Trim all alignments by reference sequence"""

    def processor(iterator):
        for refseq, seq in iterator:
            sliceobj = find_trim_slice(refseq)
            refseq = refseq[sliceobj]
            seq = seq[sliceobj]
            sliceobj = find_trim_slice(seq)
            start, end, _ = sliceobj.indices(len(seq))
            leftdots = '.' * start
            rightdots = '.' * (len(seq) - end)
            seq = seq.push_seqtext(
                leftdots + seq.seqtext[sliceobj] + rightdots,
                'boderdots()'
            )
            yield refseq, seq

    processor.command_name = 'trim-by-ref'
    processor.is_output_command = False
    return processor
