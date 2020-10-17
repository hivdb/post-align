def group_by_codons(refnas, seqnas, gap_chars='-.'):
    refcodons = []
    seqcodons = []
    bp = -1
    for refna, seqna in zip(refnas, seqnas):
        if refna not in gap_chars:
            bp = (bp + 1) % 3
            if bp == 0:
                # begin new codon
                refcodons.append([])
                seqcodons.append([])
        refcodons[-1].append(refna)
        seqcodons[-1].append(seqna)
    return refcodons, seqcodons
