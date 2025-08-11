"""Behave steps for codon alignment component tests."""

from __future__ import annotations

import os
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
import ssl
from pathlib import Path
from typing import Any

from behave import given, then, when  # type: ignore[import-untyped]


@given('sample sequences "{seqs}" and reference "{ref}"')
def given_sequences(context: Any, seqs: str, ref: str) -> None:
    """Store input sequence and reference paths.

    :param context: Behave scenario context.
    :param seqs: Path to sample sequences file.
    :param ref: Path to reference sequence file.
    """
    context.seqs_path = Path(seqs)
    context.ref_path = Path(ref)


@given('minimap2 from "{url}" is available')
def given_minimap2(context: Any, url: str) -> None:
    """Download and prepare the minimap2 binary.

    :param context: Behave scenario context.
    :param url: URL to the minimap2 tarball.
    """
    tmpdir = Path(tempfile.mkdtemp())
    archive = tmpdir / 'minimap2.tar.bz2'
    ssl_ctx = ssl.create_default_context()
    ssl_ctx.check_hostname = False
    ssl_ctx.verify_mode = ssl.CERT_NONE
    with urllib.request.urlopen(url, context=ssl_ctx) as resp:
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
    fd, path = tempfile.mkstemp(suffix='.fasta')
    os.close(fd)
    output = Path(path)
    context.output_path = output
    script = Path(__file__).with_name('run_codon_alignment.py')
    cmd = [
        sys.executable,
        str(script),
        str(context.seqs_path),
        str(context.ref_path),
        context.minimap2,
        str(output),
        opts,
        str(start),
        str(end),
    ]
    env = os.environ.copy()
    env['PATH'] = f"{Path(context.minimap2).parent}:{env['PATH']}"
    subprocess.run(cmd, check=True, env=env)


@then('the number of aligned pairs is {count:d}')
def then_check_count(context: Any, count: int) -> None:
    """Verify the expected number of sequence pairs were produced.

    :param context: Behave scenario context.
    :param count: Expected number of alignment pairs.
    """
    with open(context.output_path) as handle:
        seqs = sum(1 for line in handle if line.startswith('>')) // 2
    assert seqs == count, f"Expected {count} pairs, got {seqs}"
