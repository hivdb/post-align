"""Rust-accelerated codon alignment.

Thin Python wrapper around the ``postalign_rs`` Rust extension.  The
entire codon-alignment pipeline — boundary detection, gap gathering,
incremental IUPAC + BLOSUM62 scoring with precomputed amino-acid
tables, centre-expand search, and gap finalisation — runs in a single
Rust call (``postalign_rs.codon_align_full``).

This module exposes the same ``codon_align`` API as the pure-Python
backend in ``codon_alignment`` so the two are interchangeable.
"""

import cython  # type: ignore
import postalign_rs

from ..models import RefSeqPair, Sequence
from ..models.na_position import NAPositionList


@cython.ccall
@cython.returns(tuple)
def codon_align(
    refseq: Sequence,
    seq: Sequence,
    min_gap_distance: int,
    window_size: int,
    gap_placement_score: dict[int, dict[tuple[int, int], int]],
    ref_start: int,
    ref_end: int,
) -> RefSeqPair:
    """Perform codon-aware gap realignment via the Rust backend.

    Marshals ``NAPosition`` lists into flat arrays, delegates to
    ``postalign_rs.codon_align_full``, and reconstructs the output
    ``Sequence`` objects.

    Args:
        refseq: Reference ``Sequence`` object.
        seq: Target ``Sequence`` object.
        min_gap_distance: Merge gaps within this NA distance.
        window_size: Flanking codons used during gap search.
        gap_placement_score: ``{gap_type: {(pos, size): score}}`` map.
        ref_start: First reference position (1-based, inclusive).
        ref_end: Last reference position (1-based, inclusive).

    Returns:
        Updated ``(refseq, seq)`` pair with codon-aligned gaps.
    """
    reftext = refseq.seqtext
    seqtext = seq.seqtext

    if not isinstance(reftext, NAPositionList) or not isinstance(seqtext, NAPositionList):
        raise TypeError('codon_align requires NAPositionList seqtext')

    # Zero-copy: seqtext is already NAPositionList
    result = postalign_rs.codon_align_full_v2(
        reftext,
        seqtext,
        min_gap_distance,
        window_size,
        gap_placement_score,
        ref_start,
        ref_end,
    )

    if result is None:
        return refseq, seq

    (idx_start, idx_end, refnas, seqnas) = result

    # NAPositionList slicing + concat — no per-element Python overhead
    refseq = refseq.push_seqtext(
        reftext[:idx_start] + refnas + reftext[idx_end:],
        f'codonalign({ref_start},{ref_end})',
        0,
    )
    seq = seq.push_seqtext(
        seqtext[:idx_start] + seqnas + seqtext[idx_end:],
        f'codonalign({ref_start},{ref_end})',
        0,
    )
    return refseq, seq
