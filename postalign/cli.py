import click
from typing import (
    Union,
    TextIO,
    Optional,
    Iterable,
    List
)

from .processors import Processor
from .parsers import msa, paf, fasta, minimap2
from .models.sequence import RefSeqPair, Sequence


INPUT_FORMAT = ['MSA', 'PAF', 'MINIMAP2']


def reference_callback(
    ctx: click.Context,
    param: str,
    value: str
) -> Union[TextIO, str]:
    """Pre-process -r/--reference input"""
    retvalue: Union[TextIO, str]
    alignment_format: str = ctx.params['alignment_format']
    try:
        retvalue = open(value)
        if alignment_format == 'MSA':
            ref: Sequence = next(fasta.load(retvalue, seqtype='NA'))
            retvalue = ref.header
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception:
        if alignment_format != 'MSA':
            raise click.BadOptionUsage(
                param,
                '-r/--reference must provided as a file path '
                'if alignment is {!r}'
                .format(alignment_format),
                ctx
            )
    return retvalue


def seqs_prior_alignment_callback(
    ctx: click.Context,
    param: str,
    value: Optional[TextIO]
) -> Optional[TextIO]:
    """Pre-process -p/--seqs-prior-alignment input"""
    alignment_format: str = ctx.params['alignment_format']
    if alignment_format in ('PAF', ):
        if not value:
            raise click.BadOptionUsage(
                param,
                '-p/--seqs-prior-alignment must provided '
                'for alignment format {!r}'
                .format(alignment_format))
        return value
    elif value:
        click.echo(
            'Warning: ignore -p/--seqs-prior-alignment for '
            'alignment format {!r}'
            .format(alignment_format), err=True)
    return None


@click.group(chain=True, invoke_without_command=True)
@click.option(
    '-i', '--input-alignment',
    type=click.File('r'),
    required=True,
    help=('Input alignment file. Support formats: {}'
          .format(', '.join(INPUT_FORMAT))))
@click.option(
    '-p', '--seqs-prior-alignment',
    type=click.File('r'),
    callback=seqs_prior_alignment_callback,
    help=('FASTA sequence file prior alignment; required by '
          "'PAF' alignment format"))
@click.option(
    '-o', '--output',
    type=click.File('w'),
    required=True,
    help='Output file')
@click.option(
    '-f', '--alignment-format',
    required=True, is_eager=True,
    type=click.Choice(INPUT_FORMAT, case_sensitive=False),
    help='Input/output alignment file format')
@click.option(
    '-r', '--reference',
    type=str, required=True,
    callback=reference_callback,
    help=(
        'Header/FASTA file of the reference sequence. Will '
        'use the first sequence as reference if not specified. '
        'A file must be specified when -f/--alignment-format '
        "is not 'MSA'"
    ))
@click.option(
    '-n/-a', '--nucleotides/--amino-acids',
    default=True,
    help='The input sequences are nucleotides or amino acids')
@click.option(
    '-V/-q', '--verbose/--quiet',
    default=True,
    help='Verbose/quiet output')
@click.option(
    '--enable-profile/--disable-profile',
    default=False,
    help='Enable cProfile')
def cli(
    input_alignment: TextIO,
    seqs_prior_alignment: Optional[TextIO],
    output: TextIO,
    alignment_format: str,
    reference: Union[TextIO, str],
    nucleotides: bool,
    verbose: bool,
    enable_profile: bool
) -> None:
    pass


def call_processors(
    processors: List[Processor],
    iterator: Iterable[RefSeqPair]
) -> Iterable[str]:
    processor: Processor[Iterable[RefSeqPair]]
    for processor in processors[:-1]:
        iterator = processor(iterator)
    last_processor: Processor[Iterable[str]] = processors[-1]
    return last_processor(iterator)


def check_processors(processors: List[Processor]) -> None:
    if not processors:
        raise click.ClickException('No processor is specified')
    last_processor: Processor = processors[-1]
    if not last_processor.is_output_command:
        raise click.ClickException(
            'The last pipeline command {!r} is not an output method'
            .format(last_processor.command_name)
        )
    extra_output_commands: List[str] = []
    for processor in processors[:-1]:
        if processor.is_output_command:
            extra_output_commands.append(processor.command_name)
    if extra_output_commands:
        raise click.ClickException(
            'Following pipeline command(s) are output methods: {}'
            .format(', '.join(extra_output_commands))
        )


@cli.resultcallback()
def process_pipeline(
    processors: List[Processor],
    input_alignment: TextIO,
    seqs_prior_alignment: Optional[TextIO],
    output: TextIO,
    alignment_format: str,
    reference: Union[TextIO, str],
    nucleotides: bool,
    verbose: bool,
    enable_profile: bool
) -> None:
    seqtype: str = 'NA'
    iterator: Iterable[RefSeqPair]
    check_processors(processors)

    if not nucleotides:
        raise click.ClickException(
            'Amino acid sequences is not yet supported (--amino-acids)'
        )

    if isinstance(reference, str):
        if alignment_format == 'MSA':
            iterator = msa.load(input_alignment, reference, seqtype)
    elif alignment_format == 'PAF' and seqs_prior_alignment:
        iterator = paf.load(input_alignment, seqs_prior_alignment,
                            reference, seqtype)
    elif alignment_format == 'MINIMAP2':
        iterator = minimap2.load(input_alignment, reference, seqtype)
    else:
        raise click.ClickException(
            'Unsupport alignment format: {}'.format(alignment_format)
        )

    if enable_profile:
        import cProfile
        import pstats
        with cProfile.Profile() as profile:
            for partial in call_processors(processors, iterator):
                output.write(partial)
            ps = pstats.Stats(profile)
            ps.print_stats()
    else:
        for partial in call_processors(processors, iterator):
            output.write(partial)
