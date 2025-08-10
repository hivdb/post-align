"""Command-line interface for post-align."""

import sys
from typing import (
    Union,
    TextIO,
    Optional,
    Iterable,
    List,
    Type,
)

import typer
from rich import print

from .processor import Processor
from .parsers import msa, paf, fasta, minimap2
from .models import Message
from .models.sequence import RefSeqPair, Sequence, NAPosition


INPUT_FORMAT = ['MSA', 'PAF', 'MINIMAP2']

cli = typer.Typer(chain=True, pretty_exceptions_enable=False)


def reference_callback(
    ctx: typer.Context,
    param: typer.CallbackParam,
    value: str,
) -> Union[TextIO, str]:
    """Pre-process ``-r/--reference`` input.

    :param ctx: Typer execution context.
    :param param: Callback parameter definition.
    :param value: Provided reference value.
    :returns: Opened file handle or reference header.
    :raises typer.BadParameter: If the parameter metadata is missing.
    :raises typer.BadParameter: When file is required but missing.
    """
    if not param.name:
        raise typer.BadParameter('Internal error (reference_callback:1)')
    retvalue: Union[TextIO, str]
    alignment_format: str = ctx.params['alignment_format']
    try:
        retvalue = open(value)
        if alignment_format == 'MSA':
            ref: Sequence = next(fasta.load(retvalue, seqtype=NAPosition))
            retvalue = ref.header
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception:
        if alignment_format != 'MSA':
            raise typer.BadParameter(
                '-r/--reference must provided as a file path if alignment is '
                f"{alignment_format!r}"
            )
    return retvalue


def seqs_prior_alignment_callback(
    ctx: typer.Context,
    param: typer.CallbackParam,
    value: Optional[TextIO],
) -> Optional[TextIO]:
    """Pre-process ``-p/--seqs-prior-alignment`` input.

    :param ctx: Typer execution context.
    :param param: Callback parameter definition.
    :param value: Optional file handle.
    :returns: The provided file handle if valid, otherwise ``None``.
    :raises typer.BadParameter: If metadata is missing.
    :raises typer.BadParameter: When required input is absent.
    """
    if not param.name:
        raise typer.BadParameter(
            'Internal error (seqs_prior_alignment_callback:1)'
        )
    alignment_format: str = ctx.params['alignment_format']
    if alignment_format in ('PAF',):
        if not value:
            raise typer.BadParameter(
                '-p/--seqs-prior-alignment must provided for alignment format '
                f"{alignment_format!r}"
            )
        return value
    elif value:
        print(
            'Warning: ignore -p/--seqs-prior-alignment for alignment format '
            f"{alignment_format!r}",
            file=sys.stderr,
        )
    return None


@cli.callback(invoke_without_command=True)
def main(
    input_alignment: typer.FileText = typer.Option(
        ..., '-i', '--input-alignment',
        help=(
            'Input alignment file. Support formats: '
            f"{', '.join(INPUT_FORMAT)}"
        ),
    ),
    seqs_prior_alignment: Optional[typer.FileText] = typer.Option(
        None,
        '-p', '--seqs-prior-alignment',
        callback=seqs_prior_alignment_callback,
        help=(
            "FASTA sequence file prior alignment; required by 'PAF' "
            'alignment format'
        ),
    ),
    output: typer.FileTextWrite = typer.Option(
        ..., '-o', '--output',
        help='Output file',
    ),
    alignment_format: str = typer.Option(
        ..., '-f', '--alignment-format',
        case_sensitive=False,
        is_eager=True,
        help='Input/output alignment file format',
    ),
    reference: Union[TextIO, str] = typer.Option(
        ..., '-r', '--reference',
        callback=reference_callback,
        help=(
            'Header/FASTA file of the reference sequence. Will use the '
            'first sequence as reference if not specified. A file must be '
            "specified when -f/--alignment-format is not 'MSA'"
        ),
    ),
    nucleotides: bool = typer.Option(
        True, '-n/-a', '--nucleotides/--amino-acids',
        help='The input sequences are nucleotides or amino acids',
    ),
    verbose: bool = typer.Option(
        True, '-V/-q', '--verbose/--quiet',
        help='Verbose/quiet output',
    ),
    enable_profile: bool = typer.Option(
        False, '--enable-profile/--disable-profile',
        help='Enable cProfile',
    ),
    minimap2_opts: Optional[str] = typer.Option(
        None, '--minimap2-opts',
        help=(
            'Options to be passed to minimap2 command '
            '(when -f is MINIMAP2)'
        ),
    ),
) -> None:
    """Store common CLI options in the context."""
    pass


def call_processors(
    processors: List[Processor],
    iterator: Iterable[RefSeqPair],
    messages: List[Message]
) -> Iterable[str]:
    """Run processors sequentially.

    :param processors: Pipeline of processors.
    :param iterator: Input sequence iterator.
    :param messages: Collector for warning or error messages.
    :returns: Iterator over processed string chunks.
    """
    processor: Processor[Iterable[RefSeqPair]]
    for processor in processors[:-1]:
        iterator = processor(iterator, messages)
    last_processor: Processor[Iterable[str]] = processors[-1]
    return last_processor(iterator, messages)


def check_processors(processors: List[Processor]) -> None:
    """Validate processor pipeline configuration.

    :param processors: Sequence of processors to execute.
    :raises typer.BadParameter: If the pipeline definition is invalid.
    """
    if not processors:
        raise typer.BadParameter('No processor is specified')
    last_processor: Processor = processors[-1]
    if not last_processor.is_output_command:
        raise typer.BadParameter(
            'The last pipeline command {!r} is not an output method'
            .format(last_processor.command_name)
        )
    extra_output_commands: List[str] = []
    for processor in processors[:-1]:
        if processor.is_output_command:
            extra_output_commands.append(processor.command_name)
    if extra_output_commands:
        raise typer.BadParameter(
            'Following pipeline command(s) are output methods: {}'
            .format(', '.join(extra_output_commands))
        )


# type: ignore[attr-defined]
@cli.result_callback()  # type: ignore[attr-defined]
def process_pipeline(
    processors: List[Processor],
    input_alignment: TextIO,
    seqs_prior_alignment: Optional[TextIO],
    output: TextIO,
    alignment_format: str,
    reference: Union[TextIO, str],
    nucleotides: bool,
    verbose: bool,
    enable_profile: bool,
    minimap2_opts: Optional[str],
) -> None:
    """Execute the processor pipeline.

    :param processors: Ordered list of processors to run.
    :param input_alignment: Input alignment file handle.
    :param seqs_prior_alignment: FASTA sequence file prior to alignment.
    :param output: Output file handle.
    :param alignment_format: Alignment format identifier.
    :param reference: Reference sequence or file handle.
    :param nucleotides: ``True`` if sequences are nucleotides.
    :param verbose: Whether to output verbose messages.
    :param enable_profile: Enable profiling.
    :param minimap2_opts: Additional minimap2 options.
    :raises typer.BadParameter: On invalid configuration or unsupported format.
    """
    seqtype: Type[NAPosition] = NAPosition
    iterator: Iterable[RefSeqPair]
    check_processors(processors)
    messages: List[Message] = []

    if not nucleotides:
        raise typer.BadParameter(
            'Amino acid sequences is not yet supported (--amino-acids)'
        )

    if isinstance(reference, str):
        if alignment_format == 'MSA':
            iterator = msa.load(input_alignment, reference, seqtype)
    elif alignment_format == 'PAF' and seqs_prior_alignment:
        iterator = paf.load(input_alignment, seqs_prior_alignment,
                            reference, seqtype, messages)
    elif alignment_format == 'MINIMAP2':
        minimap2_execute = ['minimap2']
        if minimap2_opts:
            minimap2_execute.extend(minimap2_opts.split())
        iterator = minimap2.load(
            input_alignment,
            reference,
            seqtype,
            messages,
            minimap2_execute=minimap2_execute
        )
    else:
        raise typer.BadParameter(
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
