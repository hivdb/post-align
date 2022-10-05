from operator import itemgetter
from collections import defaultdict
from pafpy import PafRecord, Strand  # type: ignore
from typing import Type, Iterable, TextIO, List, Dict, Tuple, Set

from ..models import Sequence, Position, RefSeqPair
from ..utils.cigar import CIGAR

from . import fasta


def insert_unaligned_region(
    reftext: List[Position],
    seqtext: List[Position],
    orig_seqtext: List[Position],
    seqtype: Type[Position],
    align1_ref_end: int,
    align1_seq_end: int,
    align2_ref_start: int,
    align2_seq_start: int
) -> None:
    """populate unaligned region into reftext & seqtext

    For example, seq has NAs "ABC" not aligned to ref "DEFG" since
    they are too different. The three unaligned NAs are added as below:

    reftext: ...<ALIGNMENT_1>...DEF---G...<ALIGNMENT_2>...
    seqtext: ...<ALIGNMENT_1>...---ABC-...<ALIGNMENT_2>...
    """

    unaligned_ref_size: int = align2_ref_start - align1_seq_end
    unaligned_seq_size: int = align2_seq_start - align1_seq_end
    min_unaligned_size: int = min(unaligned_ref_size, unaligned_seq_size)

    if unaligned_ref_size == 0 and unaligned_seq_size == 0:
        return

    reftext[
        align1_ref_end + min_unaligned_size:
        align1_ref_end + min_unaligned_size
    ] = seqtype.init_gaps(unaligned_seq_size)

    seqtext[
        align1_ref_end + min_unaligned_size:
        align1_ref_end + min_unaligned_size
    ] = orig_seqtext[
        align1_seq_end:
        align1_seq_end + unaligned_seq_size
    ]


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
            ] = paf_lookup[str(seq.seqid)]
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
        prev_ref_start: int = len(refseq)
        prev_seq_start: int = len(seq)
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

            insert_unaligned_region(
                final_reftext,
                final_seqtext,
                seqtext,
                seqtype,
                ref_end,
                seq_end,
                prev_ref_start,
                prev_seq_start
            )

            prev_ref_start = ref_start
            prev_seq_start = seq_start

            cigar_obj: CIGAR = CIGAR(ref_start, seq_start, cigar_text)
            reftext, seqtext = cigar_obj.get_alignment(
                reftext, seqtext, seqtype)
            final_reftext[ref_start:ref_end] = reftext
            final_seqtext[ref_start:ref_end] = seqtext

        insert_unaligned_region(
            final_reftext,
            final_seqtext,
            seqtext,
            seqtype,
            0,
            0,
            prev_ref_start,
            prev_seq_start
        )

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
