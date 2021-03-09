import re
from more_itertools import chunked

CODON_TABLE = {
    'TTT': 'F',
    'TTC': 'F',
    'TTA': 'L',
    'TTG': 'L',

    'CTT': 'L',
    'CTC': 'L',
    'CTA': 'L',
    'CTG': 'L',

    'ATT': 'I',
    'ATC': 'I',
    'ATA': 'I',
    'ATG': 'M',

    'GTT': 'V',
    'GTC': 'V',
    'GTA': 'V',
    'GTG': 'V',

    'TCT': 'S',
    'TCC': 'S',
    'TCA': 'S',
    'TCG': 'S',

    'CCT': 'P',
    'CCC': 'P',
    'CCA': 'P',
    'CCG': 'P',

    'ACT': 'T',
    'ACC': 'T',
    'ACA': 'T',
    'ACG': 'T',

    'GCT': 'A',
    'GCC': 'A',
    'GCA': 'A',
    'GCG': 'A',

    'TAT': 'Y',
    'TAC': 'Y',

    'CAT': 'H',
    'CAC': 'H',
    'CAA': 'Q',
    'CAG': 'Q',

    'AAT': 'N',
    'AAC': 'N',
    'AAA': 'K',
    'AAG': 'K',

    'GAT': 'D',
    'GAC': 'D',
    'GAA': 'E',
    'GAG': 'E',

    'TGT': 'C',
    'TGC': 'C',
    'TGG': 'W',

    'CGT': 'R',
    'CGC': 'R',
    'CGA': 'R',
    'CGG': 'R',

    'AGT': 'S',
    'AGC': 'S',
    'AGA': 'R',
    'AGG': 'R',

    'GGT': 'G',
    'GGC': 'G',
    'GGA': 'G',
    'GGG': 'G',

    'TAA': '*',
    'TGA': '*',
    'TAG': '*',
}

REVERSE_CODON_TABLE = {}
for codon, aa in CODON_TABLE.items():
    REVERSE_CODON_TABLE.setdefault(aa, []).append(codon)


AMBIGUOUS_NAS = {
    'W': 'AT',
    'S': 'CG',
    'M': 'AC',
    'K': 'GT',
    'R': 'AG',
    'Y': 'CT',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ACT',
    'V': 'ACG',
    'N': 'ACGT'
}


def expand_ambiguous_na(na):
    return AMBIGUOUS_NAS.get(na, na)


def translate_codon(nas, fs_as='X', del_as='-'):
    nas = nas[:3]
    if del_as and len(nas) == 3 and all(na.is_gap() for na in nas):
        # codon is a in-frame deletion
        return del_as
    if fs_as and any(na.is_gap() for na in nas):
        # codon contains out-frame deletion
        return fs_as
    nas = ''.join(str(na) for na in nas)
    nas = nas.replace('-', 'N')
    if nas in CODON_TABLE:
        return CODON_TABLE[nas]
    aas = set()
    for na0 in AMBIGUOUS_NAS.get(nas[0], nas[0]):
        for na1 in AMBIGUOUS_NAS.get(nas[1], nas[1]):
            for na2 in AMBIGUOUS_NAS.get(nas[2], nas[2]):
                aas.add(CODON_TABLE[na0 + na1 + na2])
    CODON_TABLE[nas] = aas = ''.join(sorted(aas))
    return aas


def translate_codons(nas, fs_as='X', del_as='-'):
    all_aas = []
    for codon in chunked(nas, 3):
        aas = translate_codon(codon, fs_as, fs_as)
        all_aas.append(aas)
    return all_aas


def get_codons(aa):
    return REVERSE_CODON_TABLE[aa]


def compare_codon(base, target):
    # we assume that "base" contains only unambiguous NAs
    if re.search(r'[BDHVN]', target):
        # false if highly ambiguous NA were found
        return False
    for sna, tna in zip(base, target):
        if tna in 'ACGT':
            if tna != sna:
                break
        else:
            if sna not in AMBIGUOUS_NAS[tna]:
                break
    else:
        return True
    return False
