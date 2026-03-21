from collections.abc import Iterable
from typing import TextIO

import click

from .models import Message
from .models.sequence import NAPosition, RefSeqPair, Sequence
from .parsers import fasta, minimap2, msa, paf
from .processor import Processor

INPUT_FORMAT = ['MSA', 'PAF', 'MINIMAP2']


def reference_callback(
    ctx: click.Context, param: click.Option, value: str
) -> TextIO | str:
    """Pre-process -r/--reference input"""
    if not param.name:
        raise click.BadParameter('Internal error (reference_callback:1)', ctx, param)
    retvalue: TextIO | str
    alignment_format: str = ctx.params['alignment_format']
    try:
        retvalue = open(value)  # noqa: SIM115
        if alignment_format == 'MSA':
            ref: Sequence = next(fasta.load(retvalue, seqtype=NAPosition))
            retvalue = ref.header
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as exc:
        if alignment_format != 'MSA':
            raise click.BadOptionUsage(
                param.name,
                '-r/--reference must provided as a file path '
                f'if alignment is {alignment_format!r}',
                ctx,
            ) from exc
    return retvalue


def seqs_prior_alignment_callback(
    ctx: click.Context, param: click.Option, value: TextIO | None
) -> TextIO | None:
    """Pre-process -p/--seqs-prior-alignment input"""
    if not param.name:
        raise click.BadParameter(
            'Internal error (seqs_prior_alignment_callback:1)', ctx, param
        )
    alignment_format: str = ctx.params['alignment_format']
    if alignment_format in ('PAF',):
        if not value:
            raise click.BadOptionUsage(
                param.name,
                '-p/--seqs-prior-alignment must provided '
                f'for alignment format {alignment_format!r}',
            )
        return value
    elif value:
        click.echo(
            'Warning: ignore -p/--seqs-prior-alignment for '
            f'alignment format {alignment_format!r}',
            err=True,
        )
    return None


@click.group(chain=True, invoke_without_command=True)
@click.option(
    '-i',
    '--input-alignment',
    type=click.File('r'),
    required=True,
    help=('Input alignment file. Support formats: {}'.format(', '.join(INPUT_FORMAT))),
)
@click.option(
    '-p',
    '--seqs-prior-alignment',
    type=click.File('r'),
    callback=seqs_prior_alignment_callback,
    help=("FASTA sequence file prior alignment; required by 'PAF' alignment format"),
)
@click.option('-o', '--output', type=click.File('w'), required=True, help='Output file')
@click.option(
    '-f',
    '--alignment-format',
    required=True,
    is_eager=True,
    type=click.Choice(INPUT_FORMAT, case_sensitive=False),
    help='Input/output alignment file format',
)
@click.option(
    '-r',
    '--reference',
    type=str,
    required=True,
    callback=reference_callback,
    help=(
        'Header/FASTA file of the reference sequence. Will '
        'use the first sequence as reference if not specified. '
        'A file must be specified when -f/--alignment-format '
        "is not 'MSA'"
    ),
)
@click.option(
    '-n/-a',
    '--nucleotides/--amino-acids',
    default=True,
    help='The input sequences are nucleotides or amino acids',
)
@click.option('-V/-q', '--verbose/--quiet', default=True, help='Verbose/quiet output')
@click.option(
    '--enable-profile/--disable-profile', default=False, help='Enable cProfile'
)
@click.option(
    '--minimap2-opts',
    type=str,
    help=('Options to be passed to minimap2 command (when -f is MINIMAP2)'),
)
def cli(
    input_alignment: TextIO,
    seqs_prior_alignment: TextIO | None,
    output: TextIO,
    alignment_format: str,
    reference: TextIO | str,
    nucleotides: bool,
    verbose: bool,
    enable_profile: bool,
    minimap2_opts: str,
) -> None:
    pass


def call_processors(
    processors: list[Processor], iterator: Iterable[RefSeqPair], messages: list[Message]
) -> Iterable[str]:
    processor: Processor[Iterable[RefSeqPair]]
    for processor in processors[:-1]:
        iterator = processor(iterator, messages)
    last_processor: Processor[Iterable[str]] = processors[-1]
    return last_processor(iterator, messages)


def check_processors(processors: list[Processor]) -> None:
    if not processors:
        raise click.ClickException('No processor is specified')
    last_processor: Processor = processors[-1]
    if not last_processor.is_output_command:
        raise click.ClickException(
            f'The last pipeline command {last_processor.command_name!r} is not an output method'
        )
    extra_output_commands: list[str] = []
    for processor in processors[:-1]:
        if processor.is_output_command:
            extra_output_commands.append(processor.command_name)
    if extra_output_commands:
        raise click.ClickException(
            'Following pipeline command(s) are output methods: {}'.format(
                ', '.join(extra_output_commands)
            )
        )


@cli.result_callback()
def process_pipeline(
    processors: list[Processor],
    input_alignment: TextIO,
    seqs_prior_alignment: TextIO | None,
    output: TextIO,
    alignment_format: str,
    reference: TextIO | str,
    nucleotides: bool,
    verbose: bool,
    enable_profile: bool,
    minimap2_opts: str,
) -> None:
    seqtype: type[NAPosition] = NAPosition
    iterator: Iterable[RefSeqPair]
    check_processors(processors)
    messages: list[Message] = []

    if not nucleotides:
        raise click.ClickException(
            'Amino acid sequences is not yet supported (--amino-acids)'
        )

    if isinstance(reference, str):
        if alignment_format == 'MSA':
            iterator = msa.load(input_alignment, reference, seqtype)
    elif alignment_format == 'PAF' and seqs_prior_alignment:
        iterator = paf.load(
            input_alignment, seqs_prior_alignment, reference, seqtype, messages
        )
    elif alignment_format == 'MINIMAP2':
        minimap2_execute = ['minimap2']
        if minimap2_opts:
            minimap2_execute.extend(minimap2_opts.split())
        iterator = minimap2.load(
            input_alignment,
            reference,
            seqtype,
            messages,
            minimap2_execute=minimap2_execute,
        )
    else:
        raise click.ClickException(
            f'Unsupport alignment format: {alignment_format}'
        )

    if enable_profile:
        import cProfile
        import pstats

        with cProfile.Profile() as profile:
            for partial in call_processors(processors, iterator, messages):
                output.write(partial)
            ps = pstats.Stats(profile)
            ps.print_stats()
    else:
        for partial in call_processors(processors, iterator, messages):
            output.write(partial)
