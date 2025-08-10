"""Command-line interface for post-align."""

import sys
from enum import Enum
from typing import (
    Annotated,
    TextIO,
)
from collections.abc import Iterable

import typer
from typer import FileText, FileTextWrite
from rich import print

from .processor import Processor
from .parsers import fasta, minimap2, msa, paf
from .models import Message
from .models.sequence import NAPosition, RefSeqPair, Sequence


class AlignmentFormat(str, Enum):
    """Supported alignment file formats."""

    MSA = "MSA"
    PAF = "PAF"
    MINIMAP2 = "MINIMAP2"


cli = typer.Typer(chain=True, pretty_exceptions_enable=False)


def reference_callback(
    ctx: typer.Context,
    param: typer.CallbackParam,
    value: str,
) -> TextIO | str:
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
    retvalue: TextIO | str
    alignment_format: AlignmentFormat = ctx.params['alignment_format']
    try:
        retvalue = open(value)
        if alignment_format is AlignmentFormat.MSA:
            ref: Sequence = next(fasta.load(retvalue, seqtype=NAPosition))
            retvalue = ref.header
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception:
        if alignment_format is not AlignmentFormat.MSA:
            raise typer.BadParameter(
                '-r/--reference must provided as a file path if alignment is '
                f"{alignment_format.value!r}"
            )
    return retvalue


def seqs_prior_alignment_callback(
    ctx: typer.Context,
    param: typer.CallbackParam,
    value: TextIO | None,
) -> TextIO | None:
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
    alignment_format: AlignmentFormat = ctx.params['alignment_format']
    if alignment_format in (AlignmentFormat.PAF,):
        if not value:
            raise typer.BadParameter(
                '-p/--seqs-prior-alignment must provided for alignment format '
                f"{alignment_format.value!r}"
            )
        return value
    elif value:
        print(
            'Warning: ignore -p/--seqs-prior-alignment for alignment format '
            f"{alignment_format.value!r}",
            file=sys.stderr,
        )
    return None


@cli.callback(invoke_without_command=True)
def main(
    input_alignment: Annotated[
        FileText,
        typer.Option(
            ..., '-i', '--input-alignment',
            help=(
                'Input alignment file. Support formats: '
                f"{', '.join(fmt.value for fmt in AlignmentFormat)}"
            ),
        ),
    ],
    output: Annotated[
        FileTextWrite,
        typer.Option(..., '-o', '--output', help='Output file'),
    ],
    alignment_format: Annotated[
        AlignmentFormat,
        typer.Option(
            ..., '-f', '--alignment-format',
            case_sensitive=False,
            is_eager=True,
            help='Input/output alignment file format',
        ),
    ],
    reference: Annotated[
        TextIO | str,
        typer.Option(
            ..., '-r', '--reference',
            callback=reference_callback,
            help=(
                'Header/FASTA file of the reference sequence. Will use the '
                'first sequence as reference if not specified. A file must be '
                "specified when -f/--alignment-format is not 'MSA'"
            ),
        ),
    ],
    seqs_prior_alignment: Annotated[
        FileText | None,
        typer.Option(
            None,
            '-p', '--seqs-prior-alignment',
            callback=seqs_prior_alignment_callback,
            help=(
                "FASTA sequence file prior alignment; required by 'PAF' "
                'alignment format'
            ),
        ),
    ] = None,
    nucleotides: Annotated[
        bool,
        typer.Option(
            True, '-n/-a', '--nucleotides/--amino-acids',
            help='The input sequences are nucleotides or amino acids',
        ),
    ] = True,
    verbose: Annotated[
        bool,
        typer.Option(
            True, '-V/-q', '--verbose/--quiet',
            help='Verbose/quiet output',
        ),
    ] = True,
    enable_profile: Annotated[
        bool,
        typer.Option(
            False, '--enable-profile/--disable-profile',
            help='Enable cProfile',
        ),
    ] = False,
    minimap2_opts: Annotated[
        str | None,
        typer.Option(
            None, '--minimap2-opts',
            help=(
                'Options to be passed to minimap2 command '
                '(when -f is MINIMAP2)'
            ),
        ),
    ] = None,
) -> None:
    """Store common CLI options in the Typer context.

    :param input_alignment: Input alignment file.
    :param output: Output file handle.
    :param alignment_format: Alignment format identifier.
    :param reference: Header or FASTA file of the reference sequence.
    :param seqs_prior_alignment: FASTA sequence file prior to alignment;
        required by 'PAF' alignment format.
    :param nucleotides: ``True`` if sequences are nucleotides.
    :param verbose: Verbose or quiet output.
    :param enable_profile: Enable cProfile.
    :param minimap2_opts: Options passed to the minimap2 command.
    """


def call_processors(
    processors: list[Processor],
    iterator: Iterable[RefSeqPair],
    messages: list[Message]
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


def check_processors(processors: list[Processor]) -> None:
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
    extra_output_commands: list[str] = []
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
    processors: list[Processor],
    input_alignment: TextIO,
    seqs_prior_alignment: TextIO | None,
    output: TextIO,
    alignment_format: AlignmentFormat,
    reference: TextIO | str,
    nucleotides: bool,
    verbose: bool,
    enable_profile: bool,
    minimap2_opts: str | None,
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
    seqtype: type[NAPosition] = NAPosition
    iterator: Iterable[RefSeqPair]
    check_processors(processors)
    messages: list[Message] = []

    if not nucleotides:
        raise typer.BadParameter(
            'Amino acid sequences is not yet supported (--amino-acids)'
        )

    if isinstance(reference, str):
        if alignment_format is AlignmentFormat.MSA:
            iterator = msa.load(input_alignment, reference, seqtype)
    elif alignment_format is AlignmentFormat.PAF and seqs_prior_alignment:
        iterator = paf.load(input_alignment, seqs_prior_alignment,
                            reference, seqtype, messages)
    elif alignment_format is AlignmentFormat.MINIMAP2:
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
            f'Unsupport alignment format: {alignment_format.value}'
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
