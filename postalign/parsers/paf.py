from operator import attrgetter
from collections import defaultdict
from pafpy import PafRecord, Strand  # type: ignore
from typing import Type, Iterable, TextIO, List, Dict, Optional

from ..models import Sequence, Position, RefSeqPair
from ..utils.cigar import CIGAR

from . import fasta


def load(
    paffp: TextIO,
    seqs_prior_alignment: TextIO,
    reference: TextIO,
    seqtype: Type[Position]
) -> Iterable[RefSeqPair]:
    seq: Sequence
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
                    [],
                    modtext='error()',
                    start_offset=0
                ),
                seq.push_seqtext(
                    [],
                    modtext='error()',
                    start_offset=0
                )
            )
            continue
        if not pafs or any(paf.strand == Strand.Reverse for paf in pafs):
            # skip reverse strand alignment
            yield (
                refseq.push_seqtext(
                    [],
                    modtext='error()',
                    start_offset=0
                ),
                seq.push_seqtext(
                    [],
                    modtext='error()',
                    start_offset=0
                )
            )
            continue
        final_reftext: List[Position] = refseq.seqtext[:]
        final_seqtext: List[Position] = seqtype.init_gaps(len(final_reftext))
        prev_ref_start: Optional[int] = None
        prev_seq_start: Optional[int] = None
        ref_paf_params: List[str] = []
        seq_paf_params: List[str] = []
        for paf in sorted(pafs, key=attrgetter('tstart'), reverse=True):
            # scan PAF from end to begining
            reftext: List[Position] = refseq.seqtext
            seqtext: List[Position] = seq.seqtext

            ref_start: int = paf.tstart
            ref_end: int = paf.tend
            seq_start: int = paf.qstart
            seq_end: int = paf.qend
            cigar_text: str = paf.tags['cg'].value
            ref_paf_params.extend([str(ref_start), str(ref_end), cigar_text])
            seq_paf_params.extend([str(seq_start), str(seq_end), cigar_text])

            # multiple partial alignments; fill the gap with unaligned part
            if prev_ref_start is not None and prev_seq_start is not None:
                unaligned_size: int = min(
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

            cigar_obj: CIGAR = CIGAR(ref_start, seq_start, cigar_text)
            reftext, seqtext = cigar_obj.get_alignment(
                reftext, seqtext, seqtype)
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
