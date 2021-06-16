import json
# import math
import click
from textwrap import indent
from more_itertools import chunked

from ..cli import cli
from ..utils import group_by_gene_codons
from ..utils.codonutils import translate_codon

from .trim_by_ref import find_trim_slice


def group_gene_range_tuples(gene_range_tuples):
    tuples = []
    cur_group = []
    # use $$ as the end, it never goes into the results
    for one in gene_range_tuples + ('$$', ):
        if one.isdigit():
            cur_group.append(int(one))
        else:
            if cur_group:
                gene, *ranges = cur_group
                if len(ranges) % 2 != 0:
                    raise click.ClickException(
                        'Missing paired range values: <GENE>{}'.format(gene)
                    )
                ranges = tuple(chunked(ranges, 2))
                for refstart, refend in ranges:
                    if refstart < 1:
                        raise click.ClickException(
                            'argument <REF_START>:{} must be not less than 1'
                            .format(refstart)
                        )
                    if refend - 2 < refstart:
                        raise click.ClickException(
                            'no enough codon between arguments <REF_START>'
                            ':{} and <REF_END>:{}'
                            .format(refstart, refend)
                        )
                tuples.append((gene, ranges))
            cur_group = [one]
    return tuples


@cli.command('save-json')
@click.option(
    '--trim-by-seq/--no-trim-by-seq',
    default=True,
    help='Trim the reports by sequences (not ref sequence) or not'
)
@click.argument(
    'gene_range_tuples',
    nargs=-1
)
def save_json(trim_by_seq, gene_range_tuples):
    """Save as NucAmino style JSON reports

    <GENE_RANGE_TUPLES>:
        2n+1-tuples of
        <GENE> <REF_START1> <REF_END1> <REF_START2> <REF_END2> ...
        Use to calculate codon position within the protein/gene.
    """
    num_indent = 2
    gene_range_tuples = group_gene_range_tuples(gene_range_tuples)

    def processor(iterator):
        # TODO: MSA remap?
        yield '[\n'
        for idx, (refseq, seq) in enumerate(iterator):
            reftext = refseq.seqtext
            seqtext = seq.seqtext
            # refstart = None
            # refend = None
            if trim_by_seq:
                slicekey = find_trim_slice(seq)
                start, end, _ = slicekey.indices(len(seq))
                seqtext[:start].set_flag('trim_by_seq')
                seqtext[end:].set_flag('trim_by_seq')

                # left_reftext = reftext[:start]
                # right_reftext = reftext[end:]
                # left_reftext = left_reftext.remove_gaps()
                # right_reftext = right_reftext.remove_gaps()
                # refstart = math.ceil(len(left_reftext) / 3)
                # refend = -math.ceil(len(right_reftext) / 3)
                # if refend == 0:
                #     refend = None

            payload = {
                'Name': seq.headerdesc,
                'Modifiers': [[str(mod) for mod in mods]
                              for _, mods in seq.modifiers],
                'GeneReports': []
            }

            gene_codons = group_by_gene_codons(
                reftext, seqtext, gene_range_tuples)

            for gene, refcodons, seqcodons in gene_codons:
                codonpairs = []
                aligned_sites = []
                frameshifts = []
                for pos0, (refcd, seqcd) in enumerate(zip(refcodons,
                                                          seqcodons)):
                    if len(refcd) < 3:
                        # partial matched reference, skip
                        continue
                    # pos0 = (min(na.min_pos for na in refcd
                    #             if na.min_pos) - refstart) // 3
                    nalen = sum(not na.is_gap() for na in seqcd)
                    if 0 < nalen < 3:
                        ins_fs_len = 0
                        del_fs_len = 3 - nalen
                    else:
                        ins_fs_len = nalen % 3
                        del_fs_len = 0
                    codon_text = ''.join(str(na) for na in seqcd[:3])
                    if any(na.any_has_flag('trim_by_seq') for na in seqcd[:3]):
                        continue
                    codonpairs.append({
                        'Position': pos0 + 1,
                        'RefCodonText': ''.join(str(na) for na in refcd[:3]),
                        'CodonText': codon_text,
                        'RefAminoAcidText': translate_codon(refcd[:3]),
                        'AminoAcidText': translate_codon(seqcd[:3]),
                        'InsertedCodonsText': ''.join(
                            str(na) for na in seqcd[3:len(seqcd) - ins_fs_len]
                        ),
                        'IsInsertion': len(refcd) > 5,
                        'IsDeletion': codon_text == '---'
                    })
                    aligned_sites.append({
                        'PosAA': pos0 + 1,  # reference location
                        'PosNAs': [na.min_pos for na in seqcd],
                        'LengthNA': nalen
                    })
                    if ins_fs_len:
                        frameshifts.append({
                            'Position': pos0 + 1,
                            'GapLength': ins_fs_len,
                            'NucleicAcidsText': ''.join(
                                str(na) for na in seqcd[-ins_fs_len:]
                            ),
                            'IsInsertion': True
                        })
                    elif del_fs_len:
                        frameshifts.append({
                            'Position': pos0 + 1,
                            'GapLength': del_fs_len,
                            'IsInsertion': False
                        })
                mutations = [
                    cd for cd in codonpairs
                    if cd['RefAminoAcidText'] != cd['AminoAcidText'] or
                    cd['IsInsertion']
                ]

                if codonpairs:
                    payload['GeneReports'].append({
                        'Gene': gene,
                        'Report': {
                            'FirstAA': codonpairs[0]['Position'],
                            'LastAA': codonpairs[-1]['Position'],
                            'AlignedSites': aligned_sites,
                            'Mutations': mutations,
                            'FrameShifts': frameshifts
                        },
                        'Error': ''
                    })
                else:
                    payload['GeneReports'].append({
                        'Gene': gene,
                        'Report': {},
                        'Error': 'Sequence is not aligned'
                    })

            text = json.dumps(payload, indent=num_indent)
            if idx > 0:
                yield ',\n'
            yield indent(text, prefix=' ' * num_indent)
        yield '\n]\n'

    processor.command_name = 'save-json'
    processor.is_output_command = True
    return processor
