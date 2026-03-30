from collections.abc import Iterable
from io import StringIO
from pathlib import Path
from subprocess import PIPE, Popen, TimeoutExpired
from tempfile import TemporaryDirectory
from typing import TextIO

import click

from ..models import Message
from ..models.sequence import Position, RefSeqPair
from . import fasta, paf

DEFAULT_TIMEOUT = 300


def load(
    fastafp: TextIO,
    reference: TextIO,
    seqtype: type[Position],
    messages: list[Message],
    *,
    minimap2_execute: list[str] | None = None,
) -> Iterable[RefSeqPair]:
    dirname: str
    mm2_cmd: list[str] = ['minimap2'] if minimap2_execute is None else list(minimap2_execute)
    with TemporaryDirectory(prefix='postalign-minimap2-') as dirname:
        tempdir = Path(dirname)
        ref = next(iter(fasta.load(reference, seqtype, remove_gaps=True)))
        refpath = tempdir / 'target.fa'
        with refpath.open('w') as fp:
            fp.write(f'>{ref.headerdesc}\n{ref.seqtext_as_str}')
        seqpath = tempdir / 'query.fa'
        with seqpath.open('w') as fp:
            for seq in fasta.load(fastafp, seqtype, remove_gaps=True):
                fp.write(f'>{seq.seqid} {seq.headerdesc}\n{seq.seqtext_as_str}\n')
        proc = Popen(
            [
                *mm2_cmd,
                '-c',  # output CIGAR in PAF
                str(refpath),  # target.fa
                str(seqpath),
            ],  # query.fa
            stdout=PIPE,
            stderr=PIPE,
            encoding='utf-8',
        )
        try:
            # TODO: allow to specify timeout through input
            outs, errs = proc.communicate(timeout=DEFAULT_TIMEOUT)
        except TimeoutExpired:
            proc.kill()
            outs, errs = proc.communicate()
        if proc.returncode != 0:
            raise click.ClickException(f'Error happened during xecuting minimap2: {errs}')
        paffp = StringIO(outs)
        return paf.load(paffp, seqpath.open(), refpath.open(), seqtype, messages)
