from pafpy import PafRecord, Strand

from ..utils.cigar import CIGAR

from . import fasta


def load(paffp, seqs_prior_alignment, reference, seqtype):
    refseq = next(fasta.load(reference, seqtype, remove_gaps=True))
    seqs = fasta.load(seqs_prior_alignment, seqtype, remove_gaps=True)

    pafs = (paf.strip() for paf in paffp)
    pafs = (PafRecord.from_str(paf) for paf in pafs if paf)
    pafs = {paf.qname: paf for paf in pafs}

    for seq in seqs:
        reftext = refseq.seqtext
        try:
            paf = pafs[seq.header]
        except KeyError:
            # alignment for sequence is not found
            yield (
                refseq.push_seqtext(
                    '',
                    modtext='error()',
                    start_offset=0
                ),
                seq.push_seqtext(
                    '',
                    modtext='error()',
                    start_offset=0
                )
            )
            continue
        if paf.strand == Strand.Reverse:
            # skip reverse strand alignment
            yield (
                refseq.push_seqtext(
                    '',
                    modtext='error()',
                    start_offset=0
                ),
                seq.push_seqtext(
                    '',
                    modtext='error()',
                    start_offset=0
                )
            )
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
