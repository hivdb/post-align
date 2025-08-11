"""Behave steps for codon alignment component tests."""

from __future__ import annotations

import os
import ssl
import tarfile
import tempfile
import urllib.request
from pathlib import Path
from typing import Any, Callable

import inspect
import typer
from behave import given, then, when  # type: ignore[import-untyped]

if not hasattr(typer.Typer, "result_callback"):
    def result_callback(
        self: Any, *args: Any, **kwargs: Any
    ) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
            return func
        return decorator
    typer.Typer.result_callback = result_callback  # type: ignore[attr-defined]

if "multiple" not in inspect.signature(typer.Option).parameters:
    _orig_option = typer.Option

    def Option(*args: Any, **kwargs: Any) -> Any:  # type: ignore[override]
        kwargs.pop("multiple", None)
        return _orig_option(*args, **kwargs)
    typer.Option = Option  # type: ignore[assignment]

from postalign.cli import AlignmentFormat, process_pipeline
from postalign.processors.codon_alignment import codon_alignment
from postalign.processors.save_fasta import save_fasta


MINIMAP2_URL = (
    "https://github.com/lh3/minimap2/releases/download/"
    "v2.17/minimap2-2.17_x64-linux.tar.bz2"
)


@given('sample sequences "{seqs}" and reference "{ref}"')
def given_sequences(context: Any, seqs: str, ref: str) -> None:
    """Store input sequence and reference paths.

    :param context: Behave scenario context.
    :param seqs: Path to sample sequences file.
    :param ref: Path to reference sequence file.
    """
    context.seqs_path = Path(seqs)
    context.ref_path = Path(ref)


@given('minimap2 is available')
def given_minimap2(context: Any) -> None:
    """Download and prepare the minimap2 binary.

    :param context: Behave scenario context.
    """
    tmpdir = Path(tempfile.mkdtemp())
    archive = tmpdir / 'minimap2.tar.bz2'
    ssl_ctx = ssl.create_default_context()
    ssl_ctx.check_hostname = False
    ssl_ctx.verify_mode = ssl.CERT_NONE
    with urllib.request.urlopen(MINIMAP2_URL, context=ssl_ctx) as resp:
        with open(archive, 'wb') as out:
            out.write(resp.read())
    with tarfile.open(archive, 'r:bz2') as tar:
        tar.extractall(tmpdir)
    binary = next(tmpdir.glob('*/minimap2'))
    binary.chmod(0o755)
    context.minimap2 = str(binary)
    context.tmpdir = tmpdir


@when('they are codon-aligned with options "{opts}" from {start:d} to {end:d}')
def when_codon_aligned(context: Any, opts: str, start: int, end: int) -> None:
    """Run the codon-alignment pipeline via the CLI.

    :param context: Behave scenario context.
    :param opts: Options to pass to minimap2.
    :param start: Reference start position.
    :param end: Reference end position.
    """
    fd, path = tempfile.mkstemp(suffix='.fas')
    os.close(fd)
    output = Path(path)
    context.output_path = output
    processors = [
        codon_alignment(ref_start=start, ref_end=end),
        save_fasta(pairwise=True, modifiers=False),
    ]
    old_path = os.environ['PATH']
    os.environ['PATH'] = f"{Path(context.minimap2).parent}:{old_path}"
    try:
        with (
            open(context.seqs_path) as seqs,
            open(context.ref_path) as ref,
            open(output, 'w') as out,
        ):
            process_pipeline(
                processors,
                seqs,
                None,
                out,
                AlignmentFormat.MINIMAP2,
                ref,
                True,
                True,
                False,
                opts,
            )
    finally:
        os.environ['PATH'] = old_path


@then('the number of aligned pairs is {count:d}')
def then_check_count(context: Any, count: int) -> None:
    """Verify the expected number of sequence pairs were produced.

    :param context: Behave scenario context.
    :param count: Expected number of alignment pairs.
    """
    with open(context.output_path) as handle:
        seqs = sum(1 for line in handle if line.startswith('>')) // 2
    assert seqs == count, f"Expected {count} pairs, got {seqs}"
