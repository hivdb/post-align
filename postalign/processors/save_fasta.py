"""Processor to save sequences as FASTA."""

import typer
from typing import Annotated, Any, Iterable

from ..cli import cli
from ..models import RefSeqPair, Sequence

from ..processor import output_processor, Processor


@cli.command('save-fasta')
def save_fasta(
    preserve_order: Annotated[
        bool,
        typer.Option(
            '--preserve-order',
            help=(
                'Preserve original sequence input order / '
                'place ref sequence at first'
            ),
        ),
    ] = False,
    modifiers: Annotated[
        bool,
        typer.Option(
            '--modifiers/--no-modifiers',
            help=(
                'Include/exclude modification steps (modifiers) in '
                'sequence headers'
            ),
        ),
    ] = True,
    pairwise: Annotated[
        bool,
        typer.Option(
            '--pairwise/--msa',
            help='Save alignments in pairwise or MSA form',
        ),
    ] = False,
) -> Processor[Iterable[str]]:
    """Save prior post-alignment results as a FASTA file.

    :param preserve_order: Preserve original sequence order; place the
        reference sequence first when ``False``.
    :param modifiers: Include modification steps in sequence headers.
    :param pairwise: Save alignments in pairwise or MSA form.
    :returns: Processor yielding FASTA formatted strings.
    """

    @output_processor('save-fasta')
    def processor(
        iterator: Iterable[RefSeqPair],
        *args: Any
    ) -> Iterable[str]:
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
                yield '{}\n'.format(refseq.seqtext_as_str)

            if modifiers:
                yield '>{}\n'.format(seq.header_with_modifiers)
            else:
                yield '>{}\n'.format(seq.header)
            yield '{}\n'.format(seq.seqtext_as_str)

    return processor
