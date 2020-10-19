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
            codonpairs = []
            aligned_sites = []
            abs_seqoffset = seq.abs_seqstart
            for pos0, (refcd, seqcd) in enumerate(zip(refcodons, seqcodons)):
                codonpairs.append({
                    'Position': pos0 + 1,
                    'RefCodonText': ''.join(refcd[:3]),
                    'CodonText': ''.join(seqcd[:3]),
                    'InsertedCodonsText': ''.join(seqcd[3:]),
                    'IsInsertion': len(refcd) > 3,
                    'IsDeletion': seqcd == '---'
                })
                nalen = sum(na not in GAP_CHARS for na in seqcd)
                aligned_sites.append({
                    'PosAA': pos0 + 1,  # reference location
                    'PosNA': abs_seqoffset + 1,
                    'LengthNA': nalen
                })
                abs_seqoffset += nalen
            aligned_sites = aligned_sites[refstart:refend]
            codonpairs = codonpairs[refstart:refend]
            mutations = [cd for cd in codonpairs
                         if cd['RefCodonText'] != cd['CodonText'] and
                         not cd['IsInsertion']]
            text = json.dumps({
                'Name': seq.header,
                'Modifiers': [[str(mod) for mod in mods]
                              for _, mods in seq.modifiers],
                'Report': {
                    'FirstAA': codonpairs[0]['Position'],
                    'LastAA': codonpairs[-1]['Position'],
                    'AlignedSites': aligned_sites,
                    'Mutations': mutations
                }
            }, indent=num_indent)
            if idx > 0:
                yield ',\n'
            yield indent(text, prefix=' ' * num_indent)
        yield '\n]'

    processor.command_name = 'save-json'
    processor.is_output_command = True
    return processor
