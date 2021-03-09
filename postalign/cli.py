import click

from .parsers import msa, paf, fasta, minimap2


INPUT_FORMAT = ['MSA', 'PAF', 'MINIMAP2']


def reference_callback(ctx, param, value):
    """Pre-process -r/--reference input"""
    alignment_format = ctx.params['alignment_format']
    try:
        value = open(value)
        if alignment_format == 'MSA':
            ref = next(fasta.load(value, seqtype='NA'))
            value = ref.header
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
    return value


def seqs_prior_alignment_callback(ctx, param, value):
    """Pre-process -p/--seqs-prior-alignment input"""
    alignment_format = ctx.params['alignment_format']
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
        return


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
    type=str,
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
    input_alignment,
    seqs_prior_alignment,
    output,
    alignment_format,
    reference,
    nucleotides,
    verbose,
    enable_profile
):
    pass


def call_processors(processors, iterator):
    for processor in processors:
        iterator = processor(iterator)
    yield from iterator


@cli.resultcallback()
def process_pipeline(
    processors,
    input_alignment,
    seqs_prior_alignment,
    output,
    alignment_format,
    reference,
    nucleotides,
    verbose,
    enable_profile
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
    if alignment_format == 'MSA':
        iterator = msa.load(input_alignment, reference, seqtype)
    elif alignment_format == 'PAF':
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
