import click
from itertools import tee

from .parsers import fasta


INPUT_FORMAT = ['FASTA']


@click.group(chain=True, invoke_without_command=True)
@click.option(
    '-i', '--input-alignment',
    type=click.File('r'),
    required=True,
    help=('Input alignment file. Support formats: {}'
          .format(', '.join(INPUT_FORMAT))))
@click.option(
    '-o', '--output',
    type=click.File('w'),
    required=True,
    help='Output file')
@click.option(
    '-f', '--alignment-format',
    required=True,
    type=click.Choice(INPUT_FORMAT, case_sensitive=False),
    help='Input/output alignment file format')
@click.option(
    '-r', '--reference',
    type=str,
    help=(
        'Name of the reference sequence. Will use '
        'the first sequence as reference if not specified'
    ))
@click.option(
    '-n/-a', '--nucleotides/--amino-acids',
    default=True,
    help='The input sequences are nucleotides or amino acids')
def cli(
    input_alignment,
    output,
    alignment_format,
    reference,
    nucleotides
):
    pass


@cli.resultcallback()
def process_pipeline(
    processors,
    input_alignment,
    output,
    alignment_format,
    reference,
    nucleotides
):
    # TODO: verify if processors is ended with an output method
    if not processors:
        raise click.ClickException('No processor is specified')
    last_processor = processors[-1]
    if (
        not hasattr(last_processor, 'is_output_command') or
        not last_processor.is_output_command
    ):
        raise click.ClickException(
            'The last pipeline command {!r} is not an output method'
            .format(last_processor.command_name)
        )
    if not nucleotides:
        raise click.ClickException(
            'Amino acid sequences is not yet supported (--amino-acids)'
        )
    seqtype = 'NA'
    if alignment_format == 'FASTA':
        sequences = fasta.load(input_alignment, seqtype)
        ref_finder, sequences = tee(sequences, 2)
    else:
        raise click.ClickException(
            'Unsupport alignment format: {}'.format(alignment_format)
        )
    if reference:
        try:
            refseq = next(
                ref
                for ref in ref_finder
                if ref.header == reference
            )
        except StopIteration:
            raise click.ClickException(
                'Unable to locate reference {!r} (--reference)'
                .format(reference)
            )
    else:
        refseq = next(ref_finder)
    iterator = (
        (refseq, seq) for seq in sequences
        if seq != refseq
    )
    for processor in processors:
        iterator = processor(iterator)
    for partial in iterator:
        output.write(partial)
