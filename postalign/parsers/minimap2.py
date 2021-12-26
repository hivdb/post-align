import click
from io import StringIO
from pathlib import Path
from subprocess import Popen, TimeoutExpired, PIPE
from tempfile import TemporaryDirectory
from typing import TextIO, List, Iterable

from . import fasta, paf
from ..models.sequence import RefSeqPair

DEFAULT_TIMEOUT = 300


def load(
    fastafp: TextIO,
    reference: TextIO,
    seqtype: str,
    *,
    minimap2_execute: List[str] = ['minimap2']
) -> Iterable[RefSeqPair]:
    dirname: str
    minimap2_execute = [*minimap2_execute]  # copy in case of overwriting
    with TemporaryDirectory(prefix='postalign-minimap2-') as dirname:
        tempdir = Path(dirname)
        ref = list(fasta.load(reference, seqtype, remove_gaps=True))[0]
        refpath = tempdir / 'target.fa'
        with refpath.open('w') as fp:
            fp.write('>{}\n{}'.format(
                ref.headerdesc,
                str(ref.seqtext)
            ))
        seqpath = tempdir / 'query.fa'
        with seqpath.open('w') as fp:
            for seq in fasta.load(fastafp, seqtype, remove_gaps=True):
                fp.write('>{}\n{}\n'.format(
                    seq.headerdesc,
                    str(seq.seqtext)
                ))
        proc = Popen(
            [*minimap2_execute,
             '-c',           # output CIGAR in PAF
             str(refpath),   # target.fa
             str(seqpath)],  # query.fa
            stdout=PIPE,
            stderr=PIPE,
            encoding='utf-8'
        )
        try:
            # TODO: allow to specify timeout through input
            outs, errs = proc.communicate(timeout=DEFAULT_TIMEOUT)
        except TimeoutExpired:
            proc.kill()
            outs, errs = proc.communicate()
        if proc.returncode != 0:
            raise click.ClickException(
                'Error happened during xecuting minimap2: {}'
                .format(errs)
            )
        paffp = StringIO(outs)
        return paf.load(paffp, seqpath.open(), refpath.open(), seqtype)
