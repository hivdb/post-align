import re

from ..cli import cli


@cli.command('trim-by-ref')
def trim_by_ref():
    """Trim all alignments by reference sequence"""

    LEFT_TRIM_PATTERN = re.compile(r'^[.-]+')
    RIGHT_TRIM_PATTERN = re.compile(r'[.-]+$')

    def processor(iterator):
        trimed_refseq = None
        for refseq, seq in iterator:
            left_trim = 0
            match = LEFT_TRIM_PATTERN.search(refseq.seqtext)
            if match:
                left_trim = match.span()[1]
            right_trim = 0
            match = RIGHT_TRIM_PATTERN.search(refseq.seqtext)
            if match:
                start, end = match.span()
                right_trim = end - start
            if trimed_refseq is None:
                trimed_refseq = refseq[left_trim:-right_trim]
            trimed_seq = seq[left_trim:-right_trim]
            yield (trimed_refseq, trimed_seq)

    processor.command_name = 'trim-by-ref'
    processor.is_output_command = False
    return processor
