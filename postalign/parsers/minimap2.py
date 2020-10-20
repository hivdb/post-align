import click
from io import StringIO
from pathlib import Path
from subprocess import Popen, TimeoutExpired, PIPE
from tempfile import TemporaryDirectory

from . import fasta, paf

DEFAULT_TIMEOUT = 300


def load(fastafp, reference, seqtype, *, minimap2_execute=['minimap2']):
    minimap2_execute = [*minimap2_execute]
    with TemporaryDirectory(prefix='postalign-minimap2-') as dirname:
        tempdir = Path(dirname)
        ref = list(fasta.load(reference, seqtype, remove_gaps=True))[0]
        refpath = tempdir / 'target.fa'
        with refpath.open('w') as fp:
            fp.write('>{}\n{}'.format(ref.header, ref.seqtext))
        seqpath = tempdir / 'query.fa'
        with seqpath.open('w') as fp:
            for seq in fasta.load(fastafp, seqtype, remove_gaps=True):
                fp.write('>{}\n{}\n'.format(seq.header, seq.seqtext))
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
