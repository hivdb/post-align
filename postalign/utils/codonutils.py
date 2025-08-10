import cython  # type: ignore
from itertools import product
from more_itertools import chunked
from ..models import NAPosition

GAP_NA: int = ord(b'-')
GAP_CODON: tuple[int, int, int] = (GAP_NA,) * 3

CODON_TABLE: dict[tuple[int, ...], tuple[int, ...]] = {
    tuple(b'TTT'): tuple(b'F'),
    tuple(b'TTC'): tuple(b'F'),
    tuple(b'TTA'): tuple(b'L'),
    tuple(b'TTG'): tuple(b'L'),

    tuple(b'CTT'): tuple(b'L'),
    tuple(b'CTC'): tuple(b'L'),
    tuple(b'CTA'): tuple(b'L'),
    tuple(b'CTG'): tuple(b'L'),

    tuple(b'ATT'): tuple(b'I'),
    tuple(b'ATC'): tuple(b'I'),
    tuple(b'ATA'): tuple(b'I'),
    tuple(b'ATG'): tuple(b'M'),

    tuple(b'GTT'): tuple(b'V'),
    tuple(b'GTC'): tuple(b'V'),
    tuple(b'GTA'): tuple(b'V'),
    tuple(b'GTG'): tuple(b'V'),

    tuple(b'TCT'): tuple(b'S'),
    tuple(b'TCC'): tuple(b'S'),
    tuple(b'TCA'): tuple(b'S'),
    tuple(b'TCG'): tuple(b'S'),

    tuple(b'CCT'): tuple(b'P'),
    tuple(b'CCC'): tuple(b'P'),
    tuple(b'CCA'): tuple(b'P'),
    tuple(b'CCG'): tuple(b'P'),

    tuple(b'ACT'): tuple(b'T'),
    tuple(b'ACC'): tuple(b'T'),
    tuple(b'ACA'): tuple(b'T'),
    tuple(b'ACG'): tuple(b'T'),

    tuple(b'GCT'): tuple(b'A'),
    tuple(b'GCC'): tuple(b'A'),
    tuple(b'GCA'): tuple(b'A'),
    tuple(b'GCG'): tuple(b'A'),

    tuple(b'TAT'): tuple(b'Y'),
    tuple(b'TAC'): tuple(b'Y'),

    tuple(b'CAT'): tuple(b'H'),
    tuple(b'CAC'): tuple(b'H'),
    tuple(b'CAA'): tuple(b'Q'),
    tuple(b'CAG'): tuple(b'Q'),

    tuple(b'AAT'): tuple(b'N'),
    tuple(b'AAC'): tuple(b'N'),
    tuple(b'AAA'): tuple(b'K'),
    tuple(b'AAG'): tuple(b'K'),

    tuple(b'GAT'): tuple(b'D'),
    tuple(b'GAC'): tuple(b'D'),
    tuple(b'GAA'): tuple(b'E'),
    tuple(b'GAG'): tuple(b'E'),

    tuple(b'TGT'): tuple(b'C'),
    tuple(b'TGC'): tuple(b'C'),
    tuple(b'TGG'): tuple(b'W'),

    tuple(b'CGT'): tuple(b'R'),
    tuple(b'CGC'): tuple(b'R'),
    tuple(b'CGA'): tuple(b'R'),
    tuple(b'CGG'): tuple(b'R'),

    tuple(b'AGT'): tuple(b'S'),
    tuple(b'AGC'): tuple(b'S'),
    tuple(b'AGA'): tuple(b'R'),
    tuple(b'AGG'): tuple(b'R'),

    tuple(b'GGT'): tuple(b'G'),
    tuple(b'GGC'): tuple(b'G'),
    tuple(b'GGA'): tuple(b'G'),
    tuple(b'GGG'): tuple(b'G'),

    tuple(b'TAA'): tuple(b'*'),
    tuple(b'TGA'): tuple(b'*'),
    tuple(b'TAG'): tuple(b'*'),
}

REVERSE_CODON_TABLE: dict[int, list[bytes]] = {}
for codon, aa in CODON_TABLE.items():
    REVERSE_CODON_TABLE.setdefault(aa[0], []).append(bytes(codon))


AMBIGUOUS_NAS: dict[int, tuple[int, ...]] = {
    ord(b'W'): tuple(b'AT'),
    ord(b'S'): tuple(b'CG'),
    ord(b'M'): tuple(b'AC'),
    ord(b'K'): tuple(b'GT'),
    ord(b'R'): tuple(b'AG'),
    ord(b'Y'): tuple(b'CT'),
    ord(b'B'): tuple(b'CGT'),
    ord(b'D'): tuple(b'AGT'),
    ord(b'H'): tuple(b'ACT'),
    ord(b'V'): tuple(b'ACG'),
    ord(b'N'): tuple(b'ACGT')
}


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _translate_codon(
    nas: tuple[int, ...],
    fs_as: tuple[int],
    del_as: tuple[int],
) -> tuple[int, ...]:
    if nas in CODON_TABLE:
        return CODON_TABLE[nas]
    len_nas = len(nas)
    if del_as and len_nas == 3 and nas == GAP_CODON:
        # codon is a in-frame deletion
        return del_as
    if fs_as and (len_nas < 3 or
                  GAP_NA in nas):
        # codon contains out-frame deletion
        return fs_as

    na: int
    cand_nas: tuple[int, ...]
    aas: set[int] = set()
    unambi_nas: list[tuple[int, ...]] = []
    for na in nas:
        if na in AMBIGUOUS_NAS:
            unambi_nas.append(AMBIGUOUS_NAS[na])
        else:
            unambi_nas.append((na, ))
    for cand_nas in product(*unambi_nas):
        aas.add(CODON_TABLE[cand_nas][0])
    aas_tuple: tuple[int, ...] = tuple(sorted(aas))
    CODON_TABLE[nas] = aas_tuple
    return aas_tuple


@cython.ccall
@cython.returns(bytes)
def translate_codon(
    nas: list[NAPosition],
    fs_as: bytes = b'X',
    del_as: bytes = b'-'
) -> bytes:
    nas = nas[:3]
    nas_bytes: bytes = NAPosition.as_bytes(nas)
    aas: tuple[int, ...] = _translate_codon(
        tuple(nas_bytes),
        tuple(fs_as),
        tuple(del_as))
    return bytes(aas)


@cython.ccall
@cython.returns(list)
def translate_codons(
    nas: list[NAPosition],
    fs_as: bytes = b'X',
    del_as: bytes = b'-'
) -> list[bytes]:
    nas_bytes: bytes = NAPosition.as_bytes(nas)
    codon: list[int]
    fs_as_tuple: tuple[int, ...] = tuple(fs_as)
    del_as_tuple: tuple[int, ...] = tuple(del_as)
    all_aas: list[bytes] = []
    for codon in chunked(nas_bytes, 3):
        aas: tuple[int, ...] = _translate_codon(
            tuple(codon),
            fs_as_tuple,
            del_as_tuple)
        all_aas.append(bytes(aas))
    return all_aas


@cython.ccall
@cython.returns(list)
def get_codons(aa: int) -> list[bytes]:
    return REVERSE_CODON_TABLE[aa]


@cython.ccall
@cython.returns(cython.bint)
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
