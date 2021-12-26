import cython  # type: ignore
from typing import Optional, Tuple, List
from ..models import NAPosition


@cython.ccall
@cython.returns(tuple)
def group_by_codons(
    refnas: List[NAPosition],
    seqnas: List[NAPosition]
) -> Tuple[
    List[List[NAPosition]],
    List[List[NAPosition]]
]:
    refcodons: List[List[NAPosition]] = []
    seqcodons: List[List[NAPosition]] = []
    lastrefcodon: Optional[List[NAPosition]] = None
    lastseqcodon: Optional[List[NAPosition]] = None
    bp = -1
    for refna, seqna in zip(refnas, seqnas):
        if not refna.is_gap:
            bp = (bp + 1) % 3
            if bp == 0:
                # begin new codon
                lastrefcodon = []
                lastseqcodon = []
                refcodons.append(lastrefcodon)
                seqcodons.append(lastseqcodon)
        if lastrefcodon is not None and lastseqcodon is not None:
            lastrefcodon.append(refna)
            lastseqcodon.append(seqna)
    return refcodons, seqcodons


@cython.ccall
@cython.returns(list)
def group_by_gene_codons(
    refnas: List[NAPosition],
    seqnas: List[NAPosition],
    gene_range_tuples: List[Tuple[str, List[Tuple[int, int]]]]
) -> List[
    Tuple[
        str,
        List[List[NAPosition]],
        List[List[NAPosition]]
    ]
]:
    gene: str
    ranges: List[Tuple[int, int]]
    results: List[
        Tuple[
            str,
            List[List[NAPosition]],
            List[List[NAPosition]]
        ]
    ] = []

    for gene, ranges in gene_range_tuples:
        refcodons = []
        seqcodons = []
        for refstart, refend in ranges:
            idxstart, idxend = NAPosition.posrange2indexrange(
                refnas, refstart, refend)
            partial_refcodons, partial_seqcodons = group_by_codons(
                refnas[idxstart:idxend],
                seqnas[idxstart:idxend]
            )
            refcodons.extend(partial_refcodons)
            seqcodons.extend(partial_seqcodons)
        results.append((gene, refcodons, seqcodons))
    return results
