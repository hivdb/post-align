def group_by_codons(refnas, seqnas):
    refcodons = []
    seqcodons = []
    lastrefcodon = None
    lastseqcodon = None
    bp = -1
    for refna, seqna in zip(refnas, seqnas):
        if not refna.is_single_gap:
            bp = (bp + 1) % 3
            if bp == 0:
                # begin new codon
                lastrefcodon = []
                lastseqcodon = []
                refcodons.append(lastrefcodon)
                seqcodons.append(lastseqcodon)
        if lastrefcodon is not None:
            lastrefcodon.append(refna)
            lastseqcodon.append(seqna)
    return refcodons, seqcodons


def group_by_gene_codons(refnas, seqnas, gene_range_tuples):
    results = []

    for gene, ranges in gene_range_tuples:
        for refstart, refend in ranges:
            idxstart, idxend = refnas.posrange2indexrange(refstart, refend)
            refcodons, seqcodons = group_by_codons(
                refnas[idxstart:idxend],
                seqnas[idxstart:idxend]
            )
            results.append((gene, refcodons, seqcodons))
    return results
