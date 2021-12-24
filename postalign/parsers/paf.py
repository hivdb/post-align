from collections import defaultdict
from pafpy import PafRecord, Strand  # type: ignore
from typing import Iterable, TextIO, List, Dict

from ..models import Sequence, NAPosition, RefSeqPair
from ..utils.cigar import CIGAR

from . import fasta


def load(
    paffp: TextIO,
    seqs_prior_alignment: TextIO,
    reference: TextIO,
    seqtype: str
) -> Iterable[RefSeqPair]:
    refseq: Sequence = next(fasta.load(reference, seqtype, remove_gaps=True))
    seqs: Iterable[Sequence] = fasta.load(
        seqs_prior_alignment, seqtype, remove_gaps=True)

    pafstr_iter = (pafstr.strip() for pafstr in paffp)
    pafrec_iter = (
        PafRecord.from_str(pafstr)
        for pafstr in pafstr_iter
        if pafstr
    )
    paf_lookup: Dict[str, List[PafRecord]] = defaultdict(list)
    for pafrec in pafrec_iter:
        paf_lookup[pafrec.qname].append(pafrec)

    for seq in seqs:
        try:
            pafs: List[PafRecord] = paf_lookup[seq.header]
        except KeyError:
            # alignment for sequence is not found
            yield (
                refseq.push_seqtext(
                    NAPosition.init_empty(),
                    modtext='error()',
                    start_offset=0
                ),
                seq.push_seqtext(
                    NAPosition.init_empty(),
                    modtext='error()',
                    start_offset=0
                )
            )
            continue
        if not pafs or any(paf.strand == Strand.Reverse for paf in pafs):
            # skip reverse strand alignment
            yield (
                refseq.push_seqtext(
                    NAPosition.init_empty(),
                    modtext='error()',
                    start_offset=0
                ),
                seq.push_seqtext(
                    NAPosition.init_empty(),
                    modtext='error()',
                    start_offset=0
                )
            )
            continue
        final_reftext = refseq.seqtext[:]
        final_seqtext = NAPosition.init_gaps(len(final_reftext))
        prev_ref_start = None
        prev_seq_start = None
        ref_paf_params = []
        seq_paf_params = []
        for paf in sorted(pafs, key=lambda paf: paf.tstart, reverse=True):
            # scan PAF from end to begining
            reftext = refseq.seqtext
            seqtext = seq.seqtext

            ref_start = paf.tstart
            ref_end = paf.tend
            seq_start = paf.qstart
            seq_end = paf.qend
            cigar_text = paf.tags['cg'].value
            ref_paf_params.extend([str(ref_start), str(ref_end), cigar_text])
            seq_paf_params.extend([str(seq_start), str(seq_end), cigar_text])

            # multiple partial alignments; fill the gap with unaligned part
            if prev_seq_start is not None:
                unaligned_size = min(
                    prev_seq_start - seq_end + 1,
                    prev_ref_start - ref_end + 1
                )
                final_seqtext[
                    ref_end:
                    ref_end + unaligned_size
                ] = seqtext[
                    seq_end:
                    seq_end + unaligned_size
                ]
            prev_ref_start = ref_start
            prev_seq_start = seq_start

            cigar_obj = CIGAR(ref_start, seq_start, cigar_text)
            reftext, seqtext = cigar_obj.get_alignment(reftext, seqtext)
            final_reftext[ref_start:ref_end] = reftext
            final_seqtext[ref_start:ref_end] = seqtext
        yield (
            refseq.push_seqtext(
                final_reftext,
                modtext='paf({})'.format(','.join(ref_paf_params)),
                start_offset=ref_start
            ),
            seq.push_seqtext(
                final_seqtext,
                modtext='paf({})'.format(','.join(seq_paf_params)),
                start_offset=seq_start
            )
        )
