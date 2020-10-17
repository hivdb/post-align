import json
import math
import click
from textwrap import indent

from ..cli import cli
from ..utils import group_by_codons

from .trim_by_ref import find_trim_slice

GAP_CHARS = '-.'


@cli.command('save-json')
@click.option(
    '--trim-by-seq/--no-trim-by-seq',
    default=True,
    help='Trim the reports by sequences (not ref sequence) or not'
)
def save_json(trim_by_seq):
    """Save as NucAmino style JSON reports"""
    num_indent = 2

    def processor(iterator):
        # TODO: MSA remap?
        yield '[\n'
        for idx, (refseq, seq) in enumerate(iterator):
            reftext = refseq.seqtext
            seqtext = seq.seqtext
            refstart = None
            refend = None
            if trim_by_seq:
                slicekey = find_trim_slice(seq)
                start, end, _ = slicekey.indices(len(seq))
                left_reftext = reftext[:start]
                right_reftext = reftext[end:]
                for gap in GAP_CHARS:
                    left_reftext = left_reftext.replace(gap, '')
                    right_reftext = right_reftext.replace(gap, '')
                refstart = math.ceil(len(left_reftext) / 3)
                refend = -math.ceil(len(right_reftext) / 3)
                if refend == 0:
                    refend = None
            refcodons, seqcodons = group_by_codons(reftext, seqtext, GAP_CHARS)
            codonpairs = [
                {
                    'position': idx + 1,
                    **dict(zip(
                        ('refCodon', 'seqCodon'),
                        (''.join(cd) for cd in codons)
                    ))
                }
                for idx, codons in enumerate(zip(refcodons, seqcodons))
            ]
            codonpairs = codonpairs[refstart:refend]
            text = json.dumps({
                'name': seq.header,
                'modifiers': str(seq.modifiers),
                'codonPairs': codonpairs
            }, indent=num_indent)
            if idx > 0:
                yield ',\n'
            yield indent(text, prefix=' ' * num_indent)
        yield '\n]'

    processor.command_name = 'save-json'
    processor.is_output_command = True
    return processor
