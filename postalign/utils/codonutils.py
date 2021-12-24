from typing import Dict, List, Set, Sequence as PySequence
from more_itertools import chunked
from ..models import NAPosition

CODON_TABLE: Dict[bytes, bytes] = {
    b'TTT': b'F',
    b'TTC': b'F',
    b'TTA': b'L',
    b'TTG': b'L',

    b'CTT': b'L',
    b'CTC': b'L',
    b'CTA': b'L',
    b'CTG': b'L',

    b'ATT': b'I',
    b'ATC': b'I',
    b'ATA': b'I',
    b'ATG': b'M',

    b'GTT': b'V',
    b'GTC': b'V',
    b'GTA': b'V',
    b'GTG': b'V',

    b'TCT': b'S',
    b'TCC': b'S',
    b'TCA': b'S',
    b'TCG': b'S',

    b'CCT': b'P',
    b'CCC': b'P',
    b'CCA': b'P',
    b'CCG': b'P',

    b'ACT': b'T',
    b'ACC': b'T',
    b'ACA': b'T',
    b'ACG': b'T',

    b'GCT': b'A',
    b'GCC': b'A',
    b'GCA': b'A',
    b'GCG': b'A',

    b'TAT': b'Y',
    b'TAC': b'Y',

    b'CAT': b'H',
    b'CAC': b'H',
    b'CAA': b'Q',
    b'CAG': b'Q',

    b'AAT': b'N',
    b'AAC': b'N',
    b'AAA': b'K',
    b'AAG': b'K',

    b'GAT': b'D',
    b'GAC': b'D',
    b'GAA': b'E',
    b'GAG': b'E',

    b'TGT': b'C',
    b'TGC': b'C',
    b'TGG': b'W',

    b'CGT': b'R',
    b'CGC': b'R',
    b'CGA': b'R',
    b'CGG': b'R',

    b'AGT': b'S',
    b'AGC': b'S',
    b'AGA': b'R',
    b'AGG': b'R',

    b'GGT': b'G',
    b'GGC': b'G',
    b'GGA': b'G',
    b'GGG': b'G',

    b'TAA': b'*',
    b'TGA': b'*',
    b'TAG': b'*',
}

REVERSE_CODON_TABLE: Dict[int, List[bytes]] = {}
for codon, aa in CODON_TABLE.items():
    REVERSE_CODON_TABLE.setdefault(aa[0], []).append(codon)


AMBIGUOUS_NAS: Dict[int, bytes] = {
    ord(b'W'): b'AT',
    ord(b'S'): b'CG',
    ord(b'M'): b'AC',
    ord(b'K'): b'GT',
    ord(b'R'): b'AG',
    ord(b'Y'): b'CT',
    ord(b'B'): b'CGT',
    ord(b'D'): b'AGT',
    ord(b'H'): b'ACT',
    ord(b'V'): b'ACG',
    ord(b'N'): b'ACGT'
}


def expand_ambiguous_na(na: int) -> bytes:
    return AMBIGUOUS_NAS.get(na, bytes([na]))


def translate_codon(
    nas: PySequence[NAPosition],
    fs_as: bytes = b'X',
    del_as: bytes = b'-'
) -> bytes:
    nas = nas[:3]
    if del_as and len(nas) == 3 and \
            NAPosition.list_contains_all_gap(nas):
        # codon is a in-frame deletion
        return del_as
    if fs_as and (len(nas) < 3 or
                  NAPosition.list_contains_any_gap(nas)):
        # codon contains out-frame deletion
        return fs_as
    nas_bytes: bytes = b''.join(bytes(na) for na in nas)
    nas_bytes = nas_bytes.replace(b'-', b'N')
    if nas_bytes in CODON_TABLE:
        return CODON_TABLE[nas_bytes]

    na0: int
    na1: int
    na2: int
    aas: Set[int] = set()
    for na0 in AMBIGUOUS_NAS.get(nas_bytes[0], [nas_bytes[0]]):
        for na1 in AMBIGUOUS_NAS.get(nas_bytes[1], [nas_bytes[1]]):
            for na2 in AMBIGUOUS_NAS.get(nas_bytes[2], [nas_bytes[2]]):
                aas.add(CODON_TABLE[bytes([na0, na1, na2])][0])
    aas_bytes = bytes(sorted(aas))
    CODON_TABLE[nas_bytes] = aas_bytes
    return aas_bytes


def translate_codons(
    nas: PySequence[NAPosition],
    fs_as: bytes = b'X',
    del_as: bytes = b'-'
) -> List[bytes]:
    codon: List[NAPosition]
    all_aas: List[bytes] = []
    for codon in chunked(nas, 3):
        aas: bytes = translate_codon(codon, fs_as, fs_as)
        all_aas.append(aas)
    return all_aas


def get_codons(aa: int) -> List[bytes]:
    return REVERSE_CODON_TABLE[aa]


def compare_codon(
    base: bytes,
    target: bytes
) -> bool:
    ha: int
    sna: int
    tna: int
    # we assume that "base" contains only unambiguous NAs
    for ha in b'BDHVN':
        if ha in target:
            # false if highly ambiguous NA were found
            return False
    for sna, tna in zip(base, target):
        if tna in b'ACGT':
            if tna != sna:
                break
        else:
            if sna not in AMBIGUOUS_NAS[tna]:
                break
    else:
        return True
    return False
