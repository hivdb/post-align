"""Cython byte-array codon alignment (T3).

Same codon_align() API as codon_alignment.py but with:
  - Option D improvements (D1 precomputed AAs, D4 center-expand)
  - Flat byte-array / int-array internals instead of NAPosition objects
  - C-array lookup tables for BLOSUM62, IUPAC, codon translation
"""
import cython  # type: ignore
from typing import Tuple, List, Set, Optional, Dict
from itertools import groupby

from ..models import Sequence, RefSeqPair, NAPosition

NOGAP: int = 0b00
REFGAP: int = 0b01
SEQGAP: int = 0b10
LEFT: int = 0b00
RIGHT: int = 0b01
GAP: int = 45    # ord('-')
DOT: int = 46    # ord('.')
FS_AA: int = 88  # ord('X')
DEL_AA: int = 45 # ord('-')

# --- IUPAC bitmask table (128 entries) ---
_IUPAC_BITS: List[int] = [0] * 128
for _c, _b in [
    (b'A', 0b00001), (b'C', 0b00010), (b'G', 0b00100), (b'T', 0b01000),
    (b'W', 0b01001), (b'S', 0b00110), (b'M', 0b00011), (b'K', 0b01100),
    (b'R', 0b00101), (b'Y', 0b01010), (b'B', 0b01110), (b'D', 0b01101),
    (b'H', 0b01011), (b'V', 0b00111), (b'N', 0b01111),
]:
    _IUPAC_BITS[_c[0]] = _b
    _IUPAC_BITS[_c[0] + 32] = _b
_IUPAC_BITS[GAP] = 0b10000
_IUPAC_BITS[DOT] = 0b10000

# --- BLOSUM62 flat table (128*128) ---
_AA_ORDER: bytes = b'ACDEFGHIKLMNPQRSTVWY'
_BL = [
 4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2,
 0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2,
 -2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3,
 -1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2,
 -2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3,
 0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3,
 -2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2,
 -1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1,
 -1,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2,
 -1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1,
 -1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1,
 -2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2,
 -1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3,
 -1,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1,
 -1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2,
 1,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-2,
 0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2,
 0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1,
 -3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2,
 -2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7,
]
_BLOSUM62: List[int] = [0] * (128 * 128)
for _i, _aa_i in enumerate(_AA_ORDER):
    for _j, _aa_j in enumerate(_AA_ORDER):
        _BLOSUM62[_aa_i * 128 + _aa_j] = _BL[_i * 20 + _j]

# --- Codon translation table ---
_CT64: List[int] = [
    70,70,76,76, 83,83,83,83, 89,89,42,42, 67,67,42,87,
    76,76,76,76, 80,80,80,80, 72,72,81,81, 82,82,82,82,
    73,73,73,77, 84,84,84,84, 78,78,75,75, 83,83,82,82,
    86,86,86,86, 65,65,65,65, 68,68,69,69, 71,71,71,71,
]
_BASE_IDX: List[int] = [-1] * 128
_BASE_IDX[84] = _BASE_IDX[116] = 0  # T/t
_BASE_IDX[67] = _BASE_IDX[99] = 1   # C/c
_BASE_IDX[65] = _BASE_IDX[97] = 2   # A/a
_BASE_IDX[71] = _BASE_IDX[103] = 3  # G/g

_EXPAND: List[Optional[Tuple[int, ...]]] = [None] * 128
for _c, _idxs in [
    (b'A',(2,)),(b'C',(1,)),(b'G',(3,)),(b'T',(0,)),
    (b'W',(2,0)),(b'S',(1,3)),(b'M',(2,1)),(b'K',(3,0)),
    (b'R',(2,3)),(b'Y',(1,0)),
    (b'B',(1,3,0)),(b'D',(2,3,0)),(b'H',(2,1,0)),(b'V',(2,1,3)),
    (b'N',(2,1,3,0)),
]:
    _EXPAND[_c[0]] = _idxs
    _EXPAND[_c[0] + 32] = _idxs


# --- Scoring ---

@cython.cfunc
@cython.inline
def _iupac_score(na_a: int, na_b: int) -> float:
    if na_a == GAP and na_b == GAP:
        return 0.0
    if na_a == na_b:
        return 1.0
    bits_a: int = _IUPAC_BITS[na_a]
    bits_b: int = _IUPAC_BITS[na_b]
    u: int = bits_a | bits_b
    x: int = bits_a ^ bits_b
    uc: int = u.bit_count()
    if uc == 0:
        return 0.0
    return -(x.bit_count() / uc)


@cython.cfunc
@cython.inline
def _blosum62_score(a: bytes, b: bytes) -> float:
    total: float = 0.0
    n: int = 0
    for aa_a in a:
        for aa_b in b:
            if aa_a == DEL_AA and aa_b == DEL_AA:
                pass
            elif aa_a == DEL_AA:
                total -= 1.0
            elif aa_b == DEL_AA:
                total -= 1.0
            elif aa_a == FS_AA and aa_b == FS_AA:
                total -= 1.0
            elif aa_a == FS_AA:
                total -= 1.0
            elif aa_b == FS_AA:
                total -= 1.0
            else:
                total += _BLOSUM62[aa_a * 128 + aa_b]
            n += 1
    return total / n if n else 0.0


@cython.cfunc
@cython.inline
def _translate_codon(c0: int, c1: int, c2: int) -> bytes:
    if c0 == GAP and c1 == GAP and c2 == GAP:
        return b'-'
    if c0 == GAP or c0 == DOT or c1 == GAP or c1 == DOT \
            or c2 == GAP or c2 == DOT:
        return b'X'
    i0: int = _BASE_IDX[c0]
    i1: int = _BASE_IDX[c1]
    i2: int = _BASE_IDX[c2]
    if i0 >= 0 and i1 >= 0 and i2 >= 0:
        return bytes([_CT64[i0 * 16 + i1 * 4 + i2]])
    e0 = _EXPAND[c0]
    e1 = _EXPAND[c1]
    e2 = _EXPAND[c2]
    if e0 is None or e1 is None or e2 is None:
        return b'X'
    aas: Set[int] = set()
    for b0 in e0:
        for b1 in e1:
            for b2 in e2:
                aas.add(_CT64[b0 * 16 + b1 * 4 + b2])
    return bytes(sorted(aas))


@cython.cfunc
@cython.inline
@cython.returns(list)
def _translate_codons(nots: List[int]) -> List[bytes]:
    result: List[bytes] = []
    i: int = 0
    length: int = len(nots)
    while i + 2 < length:
        result.append(_translate_codon(nots[i], nots[i+1], nots[i+2]))
        i += 3
    if i < length:
        result.append(b'X')
    return result


@cython.cfunc
@cython.inline
def _calc_match_score(
    my_n: List[int], other_n: List[int],
    other_aas: List[bytes], base_score: float
) -> float:
    my_aas: List[bytes] = _translate_codons(my_n)
    score: float = base_score
    i: int
    for i in range(len(my_n)):
        score += _iupac_score(my_n[i], other_n[i])
    for i in range(min(len(my_aas), len(other_aas))):
        score += _blosum62_score(my_aas[i], other_aas[i])
    return score


# --- Core helpers ---

@cython.cfunc
@cython.inline
def _find_first_gap(nots: List[int]) -> int:
    i: int
    for i in range(len(nots)):
        if nots[i] == GAP or nots[i] == DOT:
            return i
    return -1


@cython.cfunc
@cython.inline
def _any_gap(nots: List[int]) -> bool:
    n: int
    for n in nots:
        if n == GAP or n == DOT:
            return True
    return False


@cython.cfunc
@cython.inline
def _all_gap(nots: List[int]) -> bool:
    n: int
    if not nots:
        return True
    for n in nots:
        if n != GAP and n != DOT:
            return False
    return True


@cython.cfunc
@cython.inline
def _count_gaps(nots: List[int]) -> int:
    c: int = 0
    n: int
    for n in nots:
        if n == GAP or n == DOT:
            c += 1
    return c


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _separate_gaps(
    nots: List[int], poss: List[int], flgs: List[int]
) -> Tuple[
    Tuple[List[int], List[int], List[int]],
    Tuple[List[int], List[int], List[int]]
]:
    nn: List[int] = []
    np_: List[int] = []
    nf: List[int] = []
    gn: List[int] = []
    gp: List[int] = []
    gf: List[int] = []
    i: int
    for i in range(len(nots)):
        if nots[i] == GAP or nots[i] == DOT:
            gn.append(nots[i])
            gp.append(poss[i])
            gf.append(flgs[i])
        else:
            nn.append(nots[i])
            np_.append(poss[i])
            nf.append(flgs[i])
    return (nn, np_, nf), (gn, gp, gf)


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _move_gaps_to_center(
    nots: List[int], poss: List[int], flgs: List[int]
) -> Tuple[List[int], List[int], List[int]]:
    (nn, np_, nf), (gn, gp, gf) = _separate_gaps(nots, poss, flgs)
    center: int = len(nn) // 2
    nn[center:center] = gn
    np_[center:center] = gp
    nf[center:center] = gf
    return nn, np_, nf


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _remove_n_gaps(
    nots: List[int], poss: List[int], flgs: List[int], n: int
) -> Tuple[List[int], List[int], List[int]]:
    rn: List[int] = []
    rp: List[int] = []
    rf: List[int] = []
    i: int
    for i in range(len(nots)):
        if (nots[i] == GAP or nots[i] == DOT) and n > 0:
            n -= 1
        else:
            rn.append(nots[i])
            rp.append(poss[i])
            rf.append(flgs[i])
    return rn, rp, rf


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _remove_redundant_gaps(
    rn: List[int], rp: List[int], rf: List[int],
    sn: List[int], sp: List[int], sf: List[int]
) -> Tuple[
    List[int], List[int], List[int],
    List[int], List[int], List[int]
]:
    n: int = min(_count_gaps(rn), _count_gaps(sn))
    if n:
        rn, rp, rf = _remove_n_gaps(rn, rp, rf, n)
        sn, sp, sf = _remove_n_gaps(sn, sp, sf, n)
    return rn, rp, rf, sn, sp, sf


@cython.cfunc
@cython.inline
@cython.returns(list)
def _find_windows_with_gap(
    ref_n: List[int], seq_n: List[int], min_gap_distance: int
) -> List[Tuple[int, int]]:
    first_gap: int = -1
    last_gap: int = -1
    windows: List[Tuple[int, int]] = []
    i: int
    for i in range(len(ref_n)):
        rg: bool = ref_n[i] == GAP or ref_n[i] == DOT
        sg: bool = seq_n[i] == GAP or seq_n[i] == DOT
        if not rg and not sg:
            continue
        if first_gap == -1:
            first_gap = last_gap = i
        elif i - last_gap > min_gap_distance:
            windows.append((first_gap, last_gap + 1))
            first_gap = last_gap = i
        else:
            last_gap = i
    if first_gap > -1:
        windows.append((first_gap, last_gap + 1))
    return windows


# --- Codon grouping (array-based) ---

Codon = Tuple[List[int], List[int], List[int]]


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _group_by_codons(
    rn: List[int], rp: List[int], rf: List[int],
    sn: List[int], sp: List[int], sf: List[int]
) -> Tuple[List[Codon], List[Codon]]:
    ref_codons: list = []
    seq_codons: list = []
    bp: int = -1
    cr_n: Optional[List[int]] = None
    cr_p: Optional[List[int]] = None
    cr_f: Optional[List[int]] = None
    cs_n: Optional[List[int]] = None
    cs_p: Optional[List[int]] = None
    cs_f: Optional[List[int]] = None
    i: int
    for i in range(len(rn)):
        if rn[i] != GAP and rn[i] != DOT:
            bp = (bp + 1) % 3
            if bp == 0:
                cr_n, cr_p, cr_f = [], [], []
                cs_n, cs_p, cs_f = [], [], []
                ref_codons.append((cr_n, cr_p, cr_f))
                seq_codons.append((cs_n, cs_p, cs_f))
        if cr_n is not None:
            cr_n.append(rn[i])
            cr_p.append(rp[i])
            cr_f.append(rf[i])
            cs_n.append(sn[i])
            cs_p.append(sp[i])
            cs_f.append(sf[i])
    return ref_codons, seq_codons


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _flatten_codons(codons: list) -> Tuple[List[int], List[int], List[int]]:
    nots: List[int] = []
    poss: List[int] = []
    flgs: List[int] = []
    for cn, cp, cf in codons:
        nots.extend(cn)
        poss.extend(cp)
        flgs.extend(cf)
    return nots, poss, flgs


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _find_codon_trim_slice(seq_codons: list) -> Tuple[int, int]:
    """Find start/stop indices excluding leading/trailing all-gap codons."""
    start: int = 0
    stop: int = len(seq_codons)
    while start < stop:
        cn, _, _ = seq_codons[start]
        if not _all_gap(cn):
            break
        start += 1
    while stop > start:
        cn, _, _ = seq_codons[stop - 1]
        if not _all_gap(cn):
            break
        stop -= 1
    return start, stop


def _codon_pair_key(item: Tuple[int, Tuple[Codon, Codon]]) -> int:
    _, (rc, sc) = item
    if _any_gap(rc[0]):
        return REFGAP
    if _any_gap(sc[0]):
        return SEQGAP
    return NOGAP


# --- extend_codons_until_gap ---

@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _extend_codons_until_gap(
    ref_cds: list, seq_cds: list, direction: int
) -> Tuple[list, list, int]:
    if direction == LEFT:
        ref_cds = list(reversed(ref_cds))
        seq_cds = list(reversed(seq_cds))
    endidx: int = len(ref_cds)
    idx: int
    for idx in range(len(ref_cds)):
        rcn, _, _ = ref_cds[idx]
        scn, _, _ = seq_cds[idx]
        if _any_gap(rcn) or _any_gap(scn):
            endidx = idx
            break
    ref_cds = ref_cds[:endidx]
    seq_cds = seq_cds[:endidx]
    if direction == LEFT:
        ref_cds.reverse()
        seq_cds.reverse()
    return ref_cds, seq_cds, len(ref_cds)


# --- move_gap_to_codon_end ---

@cython.cfunc
@cython.inline
@cython.returns(list)
def _move_gap_to_codon_end(codons: list) -> list:
    result: list = []
    for cn, cp, cf in codons:
        nn: List[int] = []
        np_: List[int] = []
        nf: List[int] = []
        gn: List[int] = []
        gp: List[int] = []
        gf: List[int] = []
        for i in range(len(cn)):
            if cn[i] == GAP or cn[i] == DOT:
                gn.append(cn[i]); gp.append(cp[i]); gf.append(cf[i])
            else:
                nn.append(cn[i]); np_.append(cp[i]); nf.append(cf[i])
        nn.extend(gn); np_.extend(gp); nf.extend(gf)
        result.append((nn, np_, nf))
    return result


# --- center_expand_positions (D4) ---

@cython.cfunc
@cython.inline
@cython.returns(list)
def _center_expand_positions(
    center: int, scanstart: int, mynas_len: int, step: int
) -> List[int]:
    positions: List[int] = []
    center = max(scanstart, ((center + step // 2) // step) * step)
    if scanstart <= center <= mynas_len:
        positions.append(center)
    offset: int = step
    while True:
        added: bool = False
        lo: int = center - offset
        hi: int = center + offset
        if lo >= scanstart:
            positions.append(lo)
            added = True
        if hi <= mynas_len and hi != lo:
            positions.append(hi)
            added = True
        if not added:
            break
        offset += step
    return positions


# --- find_best_matches (D1 + D4) ---

@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _find_best_matches(
    my_n: List[int], my_p: List[int], my_f: List[int],
    oth_n: List[int], oth_p: List[int],
    bp1_indices: Set[int],
    gap_type: int,
    gps: Dict[Tuple[int, int], int],
    is_start: bool, is_end: bool
) -> Tuple[List[int], List[int], List[int]]:
    orig_gapidx: int = _find_first_gap(my_n)
    (nn, np_, nf), (gn, gp, gf) = _separate_gaps(my_n, my_p, my_f)
    gaplen: int = len(gn)
    max_score: Optional[Tuple[float, int, int]] = None
    best_n: Optional[List[int]] = None
    best_p: Optional[List[int]] = None
    best_f: Optional[List[int]] = None
    scanstart: int = 3 if gap_type == REFGAP else 0
    nn_len: int = len(nn)

    # D1: Precompute other-side amino acids
    other_aas: List[bytes] = _translate_codons(oth_n)

    # D4: Center-expand search order
    positions: List[int] = _center_expand_positions(
        orig_gapidx, scanstart, nn_len, 3)

    idx: int
    for idx in positions:
        tn: List[int] = nn[:]
        tp: List[int] = np_[:]
        tf: List[int] = nf[:]
        tn[idx:idx] = gn
        tp[idx:idx] = gp
        tf[idx:idx] = gf
        base_score: float = float(-gaplen)
        if is_start and idx == 0:
            base_score = 0.0
        elif is_end and idx + 3 > nn_len:
            base_score = 0.0
        score_val: float = _calc_match_score(tn, oth_n, other_aas, base_score)
        napos: int
        if gap_type == REFGAP:
            napos = np_[idx - 1]
        else:
            napos = oth_p[idx]
        if (napos, gaplen) in gps:
            score_val += gps[(napos, gaplen)]
        elif (napos, 0) in gps:
            score_val += gps[(napos, 0)]

        score: Tuple[float, int, int]
        if idx in bp1_indices:
            score = (score_val + 1, 2, -idx)
        elif idx == orig_gapidx:
            score = (score_val, 1, -idx)
        else:
            score = (score_val, 0, -idx)
        if max_score is None or score > max_score:
            max_score = score
            best_n = tn
            best_p = tp
            best_f = tf
    if best_n is None:
        return nn, np_, nf
    return best_n, best_p, best_f


# --- paired_find_best_matches ---

@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _paired_find_best_matches(
    rn: List[int], rp: List[int], rf: List[int],
    sn: List[int], sp: List[int], sf: List[int],
    gap_type: int,
    gps_all: Dict[int, Dict[Tuple[int, int], int]],
    is_seq_start: bool, is_seq_end: bool
) -> Tuple[
    List[int], List[int], List[int],
    List[int], List[int], List[int]
]:
    bp1_indices: Set[int] = set()
    bp: int = 0
    i: int
    for i in range(len(rn)):
        if rn[i] == GAP or rn[i] == DOT:
            continue
        bp = (bp + 1) % 3
        if bp == 1:
            bp1_indices.add(i)

    if gap_type == REFGAP:
        rn, rp, rf = _find_best_matches(
            rn, rp, rf, sn, sp,
            bp1_indices, gap_type, gps_all[gap_type],
            False, False)
    elif gap_type == SEQGAP:
        sn, sp, sf = _find_best_matches(
            sn, sp, sf, rn, rp,
            bp1_indices, gap_type, gps_all[gap_type],
            is_seq_start, is_seq_end)
    return rn, rp, rf, sn, sp, sf


# --- gather_gaps ---

@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _gather_gaps(
    rn: List[int], rp: List[int], rf: List[int],
    sn: List[int], sp: List[int], sf: List[int],
    min_gap_distance: int
) -> Tuple[
    List[int], List[int], List[int],
    List[int], List[int], List[int]
]:
    for ws, we in reversed(_find_windows_with_gap(rn, sn, min_gap_distance)):
        wrn = rn[ws:we]; wrp = rp[ws:we]; wrf = rf[ws:we]
        wsn = sn[ws:we]; wsp = sp[ws:we]; wsf = sf[ws:we]
        wrn, wrp, wrf, wsn, wsp, wsf = _remove_redundant_gaps(
            wrn, wrp, wrf, wsn, wsp, wsf)
        wrn, wrp, wrf = _move_gaps_to_center(wrn, wrp, wrf)
        wsn, wsp, wsf = _move_gaps_to_center(wsn, wsp, wsf)
        rn[ws:we] = wrn; rp[ws:we] = wrp; rf[ws:we] = wrf
        sn[ws:we] = wsn; sp[ws:we] = wsp; sf[ws:we] = wsf
    return rn, rp, rf, sn, sp, sf


# --- adjust_gap_placement ---

@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _adjust_gap_placement(
    ref_codons: list, seq_codons: list,
    window_size: int,
    gps: Dict[int, Dict[Tuple[int, int], int]],
    is_seq_start: bool, is_seq_end: bool
) -> Tuple[list, list]:
    trim_start, trim_stop = _find_codon_trim_slice(seq_codons)

    trimmed = list(
        enumerate(zip(ref_codons, seq_codons))
    )[trim_start:trim_stop]

    for gap_type, codonpairs_iter in groupby(trimmed, _codon_pair_key):
        if gap_type == NOGAP:
            continue

        codonpairs = list(codonpairs_iter)
        start: int = codonpairs[0][0]
        end: int = codonpairs[-1][0] + 1
        refcds: list = []
        seqcds: list = []
        for _, (rc, sc) in codonpairs:
            refcds.append(rc)
            seqcds.append(sc)

        # extend left
        ext_r, ext_s, offset = _extend_codons_until_gap(
            ref_codons[max(0, start - window_size):start],
            seq_codons[max(0, start - window_size):start],
            LEFT)
        if offset:
            refcds = ext_r + refcds
            seqcds = ext_s + seqcds
            start -= offset

        # extend right
        ext_r, ext_s, offset = _extend_codons_until_gap(
            ref_codons[end:end + window_size],
            seq_codons[end:end + window_size],
            RIGHT)
        if offset:
            refcds = refcds + ext_r
            seqcds = seqcds + ext_s
            end += offset

        wrn, wrp, wrf = _flatten_codons(refcds)
        wsn, wsp, wsf = _flatten_codons(seqcds)
        wrn, wrp, wrf, wsn, wsp, wsf = _paired_find_best_matches(
            wrn, wrp, wrf, wsn, wsp, wsf,
            gap_type, gps,
            is_seq_start and start == trim_start,
            is_seq_end and end == trim_stop)
        win_rc, win_sc = _group_by_codons(wrn, wrp, wrf, wsn, wsp, wsf)
        ref_codons[start:end] = win_rc
        seq_codons[start:end] = win_sc

    return ref_codons, seq_codons


# --- realign_gaps ---

@cython.cfunc
@cython.inline
@cython.returns(tuple)
def _realign_gaps(
    rn: List[int], rp: List[int], rf: List[int],
    sn: List[int], sp: List[int], sf: List[int],
    min_gap_distance: int, window_size: int,
    gps: Dict[int, Dict[Tuple[int, int], int]],
    is_seq_start: bool, is_seq_end: bool
) -> Tuple[
    List[int], List[int], List[int],
    List[int], List[int], List[int]
]:
    rn, rp, rf, sn, sp, sf = _gather_gaps(
        rn, rp, rf, sn, sp, sf, min_gap_distance)
    ref_codons, seq_codons = _group_by_codons(rn, rp, rf, sn, sp, sf)
    ref_codons, seq_codons = _adjust_gap_placement(
        ref_codons, seq_codons, window_size, gps,
        is_seq_start, is_seq_end)
    seq_codons = _move_gap_to_codon_end(seq_codons)
    rn, rp, rf = _flatten_codons(ref_codons)
    sn, sp, sf = _flatten_codons(seq_codons)
    return rn, rp, rf, sn, sp, sf


# --- Main entry point ---

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
    refnas: List[NAPosition] = refseq.seqtext
    seqnas: List[NAPosition] = seq.seqtext

    seq_idx_start: int = 0
    seq_idx_end: int = len(seqnas)

    ref_idx_start, ref_idx_end = NAPosition.posrange2indexrange(
        refnas, ref_start, ref_end, include_boundary_gaps=True)

    idx_start: int = (
        ref_idx_start
        if ref_idx_start > seq_idx_start
        else seq_idx_start
    )

    while True:
        test_idx_start = NAPosition.min_nongap_index(refnas, idx_start)
        if test_idx_start < 0:
            break
        if (refnas[test_idx_start].pos - ref_start) % 3 == 0:
            break
        idx_start = test_idx_start + 1

    idx_end: int = (
        ref_idx_end
        if ref_idx_end < seq_idx_end
        else seq_idx_end
    )

    if idx_start == idx_end:
        return refseq, seq

    is_seq_start: bool = idx_start <= NAPosition.min_nongap_index(seqnas)
    is_seq_end: bool = idx_end > NAPosition.max_nongap_index(seqnas)

    refnas = refnas[idx_start:idx_end]
    seqnas = seqnas[idx_start:idx_end]

    if not NAPosition.any_has_gap(refnas) and \
            not NAPosition.any_has_gap(seqnas):
        return refseq, seq

    # Convert NAPosition → flat arrays
    rn: List[int] = [na.notation for na in refnas]
    rp: List[int] = [na.pos for na in refnas]
    rf: List[int] = [na.flag for na in refnas]
    sn: List[int] = [na.notation for na in seqnas]
    sp: List[int] = [na.pos for na in seqnas]
    sf: List[int] = [na.flag for na in seqnas]

    # Run array-based realign
    rn, rp, rf, sn, sp, sf = _realign_gaps(
        rn, rp, rf, sn, sp, sf,
        min_gap_distance, window_size,
        gap_placement_score,
        is_seq_start, is_seq_end)

    # Convert back to NAPosition
    new_refnas: List[NAPosition] = [
        NAPosition(n, p, f) for n, p, f in zip(rn, rp, rf)]
    new_seqnas: List[NAPosition] = [
        NAPosition(n, p, f) for n, p, f in zip(sn, sp, sf)]

    refseq = refseq.push_seqtext(
        refseq.seqtext[:idx_start] +
        new_refnas +
        refseq.seqtext[idx_end:],
        'codonalign({},{})'.format(ref_start, ref_end), 0)
    seq = seq.push_seqtext(
        seq.seqtext[:idx_start] +
        new_seqnas +
        seq.seqtext[idx_end:],
        'codonalign({},{})'.format(ref_start, ref_end), 0)
    return refseq, seq
