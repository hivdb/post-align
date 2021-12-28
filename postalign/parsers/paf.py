from operator import itemgetter
from collections import defaultdict
from pafpy import PafRecord, Strand  # type: ignore
from typing import Type, Iterable, TextIO, List, Dict, Optional, Tuple, Set

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
    ref_start: int
    ref_end: int
    seq_start: int
    seq_end: int
    cigar_text: str
    refseq: Sequence = next(fasta.load(reference, seqtype, remove_gaps=True))
    seqs: Iterable[Sequence] = fasta.load(
        seqs_prior_alignment, seqtype, remove_gaps=True)

    pafstr_iter = (pafstr.strip() for pafstr in paffp)
    pafrec_iter = (
        PafRecord.from_str(pafstr)
        for pafstr in pafstr_iter
        if pafstr
    )
    paf_lookup: Dict[
        str,
        Set[Tuple[int, int, int, int, str]]
    ] = defaultdict(set)
    for pafrec in pafrec_iter:
        if pafrec.strand == Strand.Reverse:
            continue
        paf_lookup[pafrec.qname].add((
            pafrec.tstart,
            pafrec.tend,
            pafrec.qstart,
            pafrec.qend,
            pafrec.tags['cg'].value
        ))

    for seq in seqs:
        try:
            pafs: Set[
                Tuple[int, int, int, int, str]
            ] = paf_lookup[seq.header]
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
        if not pafs:
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
        for ref_start, ref_end, seq_start, seq_end, cigar_text in \
                sorted(pafs, key=itemgetter(0), reverse=True):
            # scan PAF from end to begining
            reftext: List[Position] = refseq.seqtext
            seqtext: List[Position] = seq.seqtext

            ref_paf_params.append(
                '{},{},{}'.format(ref_start, ref_end, cigar_text))
            seq_paf_params.append(
                '{},{},{}'.format(seq_start, seq_end, cigar_text))

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
                modtext='paf({})'.format(
                    ';'.join(ref_paf_params)
                ),
                start_offset=ref_start
            ),
            seq.push_seqtext(
                final_seqtext,
                modtext='paf({})'.format(
                    ';'.join(seq_paf_params)
                ),
                start_offset=seq_start
            )
        )
