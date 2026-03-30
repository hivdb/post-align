from collections.abc import Iterable
from typing import Any

import click

from ..cli import cli
from ..models import RefSeqPair, Sequence
from ..processor import Processor, output_processor


@cli.command('save-fasta')
@click.option(
    '--preserve-order',
    is_flag=True,
    help=('Preserve original sequence input order / place ref sequence at first'),
)
@click.option(
    '--modifiers/--no-modifiers',
    default=True,
    help=('Include/exclude modification steps (modifiers) in sequence headers'),
)
@click.option('--pairwise/--msa', default=False, help='Save alignments in pairwise or MSA form')
def save_fasta(preserve_order: bool, modifiers: bool, pairwise: bool) -> Processor[Iterable[str]]:
    """Save prior post-alignment results as a FASTA file"""

    @output_processor('save-fasta')
    def processor(iterator: Iterable[RefSeqPair], *args: Any) -> Iterable[str]:
        # TODO: MSA remap?
        idx: int
        refseq: Sequence
        seq: Sequence
        for idx, (refseq, seq) in enumerate(iterator):
            if pairwise or (not preserve_order and idx == 0) or (preserve_order and refseq.seqid + 1 == seq.seqid):
                if modifiers:
                    yield f'>{refseq.header_with_modifiers}\n'
                else:
                    yield f'>{refseq.header}\n'
                yield f'{refseq.seqtext_as_str}\n'

            if modifiers:
                yield f'>{seq.header_with_modifiers}\n'
            else:
                yield f'>{seq.header}\n'
            yield f'{seq.seqtext_as_str}\n'

    return processor
