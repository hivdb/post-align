import click
from pafpy import PafRecord, Strand

from ..utils.cigar import CIGAR

from . import fasta


def load(paffp, seqs_prior_alignment, reference, seqtype):
    refseq = next(fasta.load(reference, seqtype, remove_gaps=True))
    seqs = fasta.load(seqs_prior_alignment, seqtype, remove_gaps=True)
    seqs = {seq.header: seq for seq in seqs}

    pafs = (paf.strip() for paf in paffp)
    pafs = (PafRecord.from_str(paf) for paf in pafs if paf)

    for paf in pafs:
        reftext = refseq.seqtext
        try:
            seq = seqs[paf.qname]
        except KeyError:
            click.echo(
                'Warning: sequence for alignment {!r} '
                'is not found'
                .format(paf.qname),
                err=True)
            continue
        if paf.strand == Strand.Reverse:
            click.echo(
                'Warning: skip reverse strand alignment {!r}'
                .format(paf.qname),
                err=True)
            continue
        seqtext = seq.seqtext
        seq_start = paf.qstart
        ref_start = paf.tstart
        cigar_text = paf.tags['cg'].value
        cigar_obj = CIGAR(ref_start, seq_start, cigar_text)
        reftext, seqtext = cigar_obj.get_alignment(reftext, seqtext)
        yield (
            refseq.push_seqtext(
                reftext,
                modtext='paf({},{})'.format(ref_start, cigar_text),
                start_offset=ref_start
            ),
            seq.push_seqtext(
                seqtext,
                modtext='paf({},{})'.format(seq_start, cigar_text),
                start_offset=seq_start
            )
        )
