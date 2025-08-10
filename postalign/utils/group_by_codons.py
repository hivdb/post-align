import cython  # type: ignore
from ..models import NAPosition


@cython.ccall
@cython.returns(tuple)
def group_by_codons(
    refnas: list[NAPosition],
    seqnas: list[NAPosition]
) -> tuple[
    list[list[NAPosition]],
    list[list[NAPosition]]
]:
    refcodons: list[list[NAPosition]] = []
    seqcodons: list[list[NAPosition]] = []
    lastrefcodon: list[NAPosition] | None = None
    lastseqcodon: list[NAPosition] | None = None
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
@cython.returns(slice)
def find_codon_trim_slice(
    codons: list[list[NAPosition]]
) -> slice:
    left_trim: int = 0
    for idx, codon in enumerate(codons):
        if NAPosition.all_have_gap(codon):
            left_trim = idx + 1
        else:
            break
    codons_len: int = len(codons)
    right_trim: int = codons_len
    for idx, codon in enumerate(reversed(codons)):
        if NAPosition.all_have_gap(codon):
            right_trim = codons_len - 1 - idx
        else:
            break
    return slice(left_trim, right_trim)


@cython.ccall
@cython.returns(list)
def group_by_gene_codons(
    refnas: list[NAPosition],
    seqnas: list[NAPosition],
    gene_range_tuples: list[tuple[str, list[tuple[int, int]]]]
) -> list[
    tuple[
        str,
        list[list[NAPosition]],
        list[list[NAPosition]]
    ]
]:
    gene: str
    ranges: list[tuple[int, int]]
    results: list[
        tuple[
            str,
            list[list[NAPosition]],
            list[list[NAPosition]]
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
