"""Rust-accelerated codon alignment (T5).

Full-Rust pipeline: boundary detection, incremental scoring engine,
precomputed IUPAC table, InlineAA codon translation, and pre-allocated
work buffers — all executed in a single Rust call via codon_align_full.
"""
import cython  # type: ignore
from typing import Tuple, Dict

from ..models import Sequence, RefSeqPair, NAPosition

import postalign_rs  # Rust extension


# -----------------------------------------------------------------------
# Main entry point — same API as original
# -----------------------------------------------------------------------

@cython.ccall
@cython.returns(tuple)
def codon_align(
    refseq: Sequence,
    seq: Sequence,
    min_gap_distance: int,
    window_size: int,
    gap_placement_score: Dict[int, Dict[Tuple[int, int], int]],
    ref_start: int,
    ref_end: int
) -> RefSeqPair:
    # Extract flat arrays from NAPosition lists (one pass each)
    reftext = refseq.seqtext
    seqtext = seq.seqtext

    ref_notations = [na.notation for na in reftext]
    ref_positions = [na.pos for na in reftext]
    ref_flags = [na.flag for na in reftext]
    seq_notations = [na.notation for na in seqtext]
    seq_positions = [na.pos for na in seqtext]
    seq_flags = [na.flag for na in seqtext]

    # Single Rust call: boundary detection + realign_gaps
    result = postalign_rs.codon_align_full(
        ref_notations, ref_positions, ref_flags,
        seq_notations, seq_positions, seq_flags,
        min_gap_distance,
        window_size,
        gap_placement_score,
        ref_start,
        ref_end,
    )

    if result is None:
        # No alignment needed (no gaps or empty window)
        return refseq, seq

    (idx_start, idx_end,
     out_ref_n, out_ref_p, out_ref_f,
     out_seq_n, out_seq_p, out_seq_f) = result

    # Reconstruct NAPosition lists from Rust output
    refnas = [
        NAPosition(n, p, f)
        for n, p, f in zip(out_ref_n, out_ref_p, out_ref_f)
    ]
    seqnas = [
        NAPosition(n, p, f)
        for n, p, f in zip(out_seq_n, out_seq_p, out_seq_f)
    ]

    # Save "codon aligned" refseq and seq
    refseq = refseq.push_seqtext(
        reftext[:idx_start] +
        refnas +
        reftext[idx_end:],
        'codonalign({},{})'.format(ref_start, ref_end), 0)
    seq = seq.push_seqtext(
        seqtext[:idx_start] +
        seqnas +
        seqtext[idx_end:],
        'codonalign({},{})'.format(ref_start, ref_end), 0)
    return refseq, seq
