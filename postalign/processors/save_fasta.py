import click
from typing import Iterable

from ..cli import cli
from ..models import RefSeqPair, Sequence

from ..processor import output_processor, Processor


@cli.command('save-fasta')
@click.option(
    '--preserve-order',
    is_flag=True,
    help=(
        'Preserve original sequence input order / '
        'place ref sequence at first'))
@click.option(
    '--modifiers/--no-modifiers',
    default=True,
    help=(
        'Include/exclude modification steps (modifiers) '
        'in sequence headers'))
@click.option(
    '--pairwise/--msa',
    default=False,
    help='Save alignments in pairwise or MSA form')
def save_fasta(
    preserve_order: bool,
    modifiers: bool,
    pairwise: bool
) -> Processor[Iterable[str]]:
    """Save prior post-alignment results as a FASTA file"""

    @output_processor('save-fasta')
    def processor(iterator: Iterable[RefSeqPair]) -> Iterable[str]:
        # TODO: MSA remap?
        idx: int
        refseq: Sequence
        seq: Sequence
        for idx, (refseq, seq) in enumerate(iterator):
            if (
                pairwise or
                (not preserve_order and idx == 0) or
                (preserve_order and refseq.seqid + 1 == seq.seqid)
            ):
                if modifiers:
                    yield '>{}\n'.format(refseq.header_with_modifiers)
                else:
                    yield '>{}\n'.format(refseq.header)
                yield '{}\n'.format(refseq.seqtext)

            if modifiers:
                yield '>{}\n'.format(seq.header_with_modifiers)
            else:
                yield '>{}\n'.format(seq.header)
            yield '{}\n'.format(seq.seqtext)

    return processor
