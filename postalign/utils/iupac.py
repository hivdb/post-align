import cython  # type: ignore
from typing import Dict, Set

IUPAC: Dict[int, Set[int]] = {
    ord(b'A'): set(b'A'),
    ord(b'C'): set(b'C'),
    ord(b'G'): set(b'G'),
    ord(b'T'): set(b'T'),
    ord(b'W'): set(b'AT'),
    ord(b'S'): set(b'CG'),
    ord(b'M'): set(b'AC'),
    ord(b'K'): set(b'GT'),
    ord(b'R'): set(b'AG'),
    ord(b'Y'): set(b'CT'),
    ord(b'B'): set(b'CGT'),
    ord(b'D'): set(b'AGT'),
    ord(b'H'): set(b'ACT'),
    ord(b'V'): set(b'ACG'),
    ord(b'N'): set(b'ACGT'),
    ord(b'-'): set(b'-')
}


@cython.ccall
def iupac_score(
    na_a: int,
    na_b: int,
    del_as: int = ord(b'-')
) -> float:
    if na_a == na_b == del_as:
        # both are in-frame deletions, no penalty applied
        return 0
    elif na_a == na_b:
        return 1
    else:
        expand_a: Set[int] = IUPAC[na_a]
        expand_b: Set[int] = IUPAC[na_b]
        return - len(expand_a ^ expand_b) / len(expand_a | expand_b)
