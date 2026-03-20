//! Scoring tables and functions: BLOSUM62, IUPAC, codon translation.
//!
//! T5 optimizations: precomputed IUPAC lookup table, inline codon
//! translation with [u8;4] (no heap), precomputed BLOSUM62 table.

// ---------------------------------------------------------------------------
// IUPAC nucleotide scoring — precomputed 128×128 table
// ---------------------------------------------------------------------------

/// Return bitmask for IUPAC nucleotide: A=1, C=2, G=4, T=8, gap=16
const fn iupac_bits(na: u8) -> u32 {
    match na {
        b'A' | b'a' => 0b00001,
        b'C' | b'c' => 0b00010,
        b'G' | b'g' => 0b00100,
        b'T' | b't' => 0b01000,
        b'W' | b'w' => 0b01001, // AT
        b'S' | b's' => 0b00110, // CG
        b'M' | b'm' => 0b00011, // AC
        b'K' | b'k' => 0b01100, // GT
        b'R' | b'r' => 0b00101, // AG
        b'Y' | b'y' => 0b01010, // CT
        b'B' | b'b' => 0b01110, // CGT
        b'D' | b'd' => 0b01101, // AGT
        b'H' | b'h' => 0b01011, // ACT
        b'V' | b'v' => 0b00111, // ACG
        b'N' | b'n' => 0b01111, // ACGT
        b'-' | b'.' => 0b10000, // gap
        _ => 0,
    }
}

/// Build precomputed IUPAC score table at compile time.
/// Indexed by [na_a * 128 + na_b]. Stores score × 1000 as i32
/// to allow const fn (no f64 ops in const). We convert at lookup.
const fn build_iupac_table() -> [i32; 128 * 128] {
    let mut table = [0i32; 128 * 128];
    // All IUPAC characters we care about
    let chars: &[u8] = b"ACGTWSMKRYBDHVNacgtwsmkrybdhvn-.";
    let mut i = 0;
    while i < chars.len() {
        let mut j = 0;
        while j < chars.len() {
            let a = chars[i];
            let b = chars[j];
            let idx = a as usize * 128 + b as usize;
            if a == b'-' && b == b'-' {
                table[idx] = 0;
            } else if a == b {
                table[idx] = 1000;
            } else {
                let bits_a = iupac_bits(a);
                let bits_b = iupac_bits(b);
                let union = bits_a | bits_b;
                let xor = bits_a ^ bits_b;
                let union_count = union.count_ones();
                let xor_count = xor.count_ones();
                if union_count == 0 {
                    table[idx] = 0;
                } else {
                    // -(xor/union) stored as *1000
                    table[idx] = -((xor_count * 1000 / union_count) as i32);
                }
            }
            j += 1;
        }
        i += 1;
    }
    table
}

static IUPAC_TABLE: [i32; 128 * 128] = build_iupac_table();

#[inline(always)]
pub fn iupac_score(na_a: u8, na_b: u8) -> f64 {
    IUPAC_TABLE[na_a as usize * 128 + na_b as usize] as f64 / 1000.0
}

// ---------------------------------------------------------------------------
// BLOSUM62
// ---------------------------------------------------------------------------

/// Amino acid order: ACDEFGHIKLMNPQRSTVWY (20 standard)
const AA_ORDER: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

#[rustfmt::skip]
const BLOSUM62_FLAT: [i8; 400] = [
//   A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
     4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2, // A
     0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2, // C
    -2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3, // D
    -1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2, // E
    -2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3, // F
     0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3, // G
    -2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2, // H
    -1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1, // I
    -1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2, // K
    -1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1, // L
    -1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1, // M
    -2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2, // N
    -1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3, // P
    -1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1, // Q
    -1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2, // R
     1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2, // S
     0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2, // T
     0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1, // V
    -3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2, // W
    -2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7, // Y
];

/// Build 128×128 lookup table from the flat 20×20 matrix.
const fn build_blosum62_table() -> [i8; 128 * 128] {
    let mut table = [0i8; 128 * 128];
    let aa = AA_ORDER;
    let mut i = 0;
    while i < 20 {
        let mut j = 0;
        while j < 20 {
            table[aa[i] as usize * 128 + aa[j] as usize] = BLOSUM62_FLAT[i * 20 + j];
            j += 1;
        }
        i += 1;
    }
    table
}

static BLOSUM62_TABLE: [i8; 128 * 128] = build_blosum62_table();

#[inline]
fn blosum62_lookup(aa_a: u8, aa_b: u8) -> i8 {
    BLOSUM62_TABLE[aa_a as usize * 128 + aa_b as usize]
}

/// Score two amino-acid byte slices (may be multi-AA for ambiguous codons).
/// Matches Python `blosum62_score()` exactly.
pub fn blosum62_score(a: &[u8], b: &[u8]) -> f64 {
    const FS: u8 = b'X';
    const DEL: u8 = b'-';
    const FS_PEN: f64 = -1.0;
    const DEL_PEN: f64 = -1.0;

    let mut total_score: f64 = 0.0;
    let mut total_n: usize = 0;

    for &aa_a in a {
        for &aa_b in b {
            if aa_a == DEL && aa_b == DEL {
                // both in-frame deletions, no penalty
            } else if aa_a == DEL {
                total_score += DEL_PEN;
            } else if aa_b == DEL {
                total_score += DEL_PEN;
            } else if aa_a == FS && aa_b == FS {
                total_score += FS_PEN; // (a_fs + b_fs) / 2 = -1
            } else if aa_a == FS {
                total_score += FS_PEN;
            } else if aa_b == FS {
                total_score += FS_PEN;
            } else {
                total_score += blosum62_lookup(aa_a, aa_b) as f64;
            }
            total_n += 1;
        }
    }
    if total_n == 0 {
        0.0
    } else {
        total_score / total_n as f64
    }
}

// ---------------------------------------------------------------------------
// Codon translation
// ---------------------------------------------------------------------------

/// Standard genetic code, indexed by (base1*16 + base2*4 + base3)
/// where T=0, C=1, A=2, G=3.
#[rustfmt::skip]
const CODON_TABLE_64: [u8; 64] = [
    b'F', b'F', b'L', b'L',  b'S', b'S', b'S', b'S', // TT*, TC*
    b'Y', b'Y', b'*', b'*',  b'C', b'C', b'*', b'W', // TA*, TG*
    b'L', b'L', b'L', b'L',  b'P', b'P', b'P', b'P', // CT*, CC*
    b'H', b'H', b'Q', b'Q',  b'R', b'R', b'R', b'R', // CA*, CG*
    b'I', b'I', b'I', b'M',  b'T', b'T', b'T', b'T', // AT*, AC*
    b'N', b'N', b'K', b'K',  b'S', b'S', b'R', b'R', // AA*, AG*
    b'V', b'V', b'V', b'V',  b'A', b'A', b'A', b'A', // GT*, GC*
    b'D', b'D', b'E', b'E',  b'G', b'G', b'G', b'G', // GA*, GG*
];

#[inline]
fn base_to_idx(b: u8) -> Option<usize> {
    match b {
        b'T' | b't' => Some(0),
        b'C' | b'c' => Some(1),
        b'A' | b'a' => Some(2),
        b'G' | b'g' => Some(3),
        _ => None,
    }
}

/// Expand an IUPAC base to its unambiguous components.
fn expand_base(b: u8) -> &'static [u8] {
    match b {
        b'A' | b'a' => &[b'A'],
        b'C' | b'c' => &[b'C'],
        b'G' | b'g' => &[b'G'],
        b'T' | b't' => &[b'T'],
        b'W' | b'w' => &[b'A', b'T'],
        b'S' | b's' => &[b'C', b'G'],
        b'M' | b'm' => &[b'A', b'C'],
        b'K' | b'k' => &[b'G', b'T'],
        b'R' | b'r' => &[b'A', b'G'],
        b'Y' | b'y' => &[b'C', b'T'],
        b'B' | b'b' => &[b'C', b'G', b'T'],
        b'D' | b'd' => &[b'A', b'G', b'T'],
        b'H' | b'h' => &[b'A', b'C', b'T'],
        b'V' | b'v' => &[b'A', b'C', b'G'],
        b'N' | b'n' => &[b'A', b'C', b'G', b'T'],
        _ => &[],
    }
}

const GAP: u8 = b'-';
const DOT: u8 = b'.';
const FS_AA: u8 = b'X';
const DEL_AA: u8 = b'-';

// ---------------------------------------------------------------------------
// Inline codon translation — zero heap allocation
// ---------------------------------------------------------------------------

/// Inline amino acid result: up to 4 AAs stored on the stack.
#[derive(Clone, Copy)]
pub struct InlineAA {
    pub data: [u8; 4],
    pub len: u8,
}

impl InlineAA {
    #[inline(always)]
    pub fn as_slice(&self) -> &[u8] {
        &self.data[..self.len as usize]
    }
}

/// Translate a 3-byte codon to amino acid(s) without heap allocation.
#[inline]
pub fn translate_codon_inline(codon: &[u8]) -> InlineAA {
    let len = codon.len();
    // All-gap codon → in-frame deletion
    if len == 3 && codon[0] == GAP && codon[1] == GAP && codon[2] == GAP {
        return InlineAA { data: [DEL_AA, 0, 0, 0], len: 1 };
    }
    // Frameshift: <3 bases or contains a gap
    if len < 3 || codon[0] == GAP || codon[0] == DOT
              || codon[1] == GAP || codon[1] == DOT
              || codon[2] == GAP || codon[2] == DOT {
        return InlineAA { data: [FS_AA, 0, 0, 0], len: 1 };
    }
    // Standard (unambiguous) codon
    if let (Some(i0), Some(i1), Some(i2)) =
        (base_to_idx(codon[0]), base_to_idx(codon[1]), base_to_idx(codon[2]))
    {
        let aa = CODON_TABLE_64[i0 * 16 + i1 * 4 + i2];
        return InlineAA { data: [aa, 0, 0, 0], len: 1 };
    }
    // Ambiguous codon
    let e0 = expand_base(codon[0]);
    let e1 = expand_base(codon[1]);
    let e2 = expand_base(codon[2]);
    let mut result = InlineAA { data: [0; 4], len: 0 };
    for &b0 in e0 {
        for &b1 in e1 {
            for &b2 in e2 {
                if let (Some(i0), Some(i1), Some(i2)) =
                    (base_to_idx(b0), base_to_idx(b1), base_to_idx(b2))
                {
                    let aa = CODON_TABLE_64[i0 * 16 + i1 * 4 + i2];
                    let n = result.len as usize;
                    let mut found = false;
                    let mut k = 0;
                    while k < n {
                        if result.data[k] == aa { found = true; break; }
                        k += 1;
                    }
                    if !found && n < 4 {
                        result.data[n] = aa;
                        result.len += 1;
                    }
                }
            }
        }
    }
    // Sort the result for consistency
    let n = result.len as usize;
    if n > 1 {
        let s = &mut result.data[..n];
        s.sort_unstable();
    }
    result
}

/// Vec-based translate_codon for backward compatibility with T4.
pub fn translate_codon(codon: &[u8]) -> Vec<u8> {
    let inline = translate_codon_inline(codon);
    inline.as_slice().to_vec()
}

/// Translate a flat nucleotide array into per-codon amino acid lists.
pub fn translate_codons(nas: &[u8]) -> Vec<Vec<u8>> {
    nas.chunks(3).map(|chunk| translate_codon(chunk)).collect()
}

/// Translate codons into inline AAs (no heap per codon).
pub fn translate_codons_inline(nas: &[u8]) -> Vec<InlineAA> {
    nas.chunks(3).map(|chunk| translate_codon_inline(chunk)).collect()
}

// ---------------------------------------------------------------------------
// Combined match score (IUPAC + BLOSUM62)
// ---------------------------------------------------------------------------

use crate::align::NaPos;

/// Compute match score with precomputed other-side amino acids (D1).
/// Kept for backward compatibility with T4 align.rs.
pub fn calc_match_score_precomputed(
    mynas: &[NaPos],
    othernas: &[NaPos],
    other_aas: &[Vec<u8>],
    base_score: f64,
) -> f64 {
    let my_bytes: Vec<u8> = mynas.iter().map(|n| n.notation).collect();
    let my_aas = translate_codons(&my_bytes);
    let mut score = base_score;

    // IUPAC positional score
    for (myna, otherna) in mynas.iter().zip(othernas.iter()) {
        score += iupac_score(myna.notation, otherna.notation);
    }

    // BLOSUM62 codon-level score
    for (myaa, otheraa) in my_aas.iter().zip(other_aas.iter()) {
        score += blosum62_score(myaa, otheraa);
    }

    score
}

// ---------------------------------------------------------------------------
// T5: InlineAA-based BLOSUM62 scoring (no heap)
// ---------------------------------------------------------------------------

/// Score two InlineAA codon translations.
#[inline]
pub fn blosum62_score_inline(a: &InlineAA, b: &InlineAA) -> f64 {
    const FS: u8 = b'X';
    const DEL: u8 = b'-';
    const FS_PEN: f64 = -1.0;
    const DEL_PEN: f64 = -1.0;

    let mut total_score: f64 = 0.0;
    let mut total_n: usize = 0;
    let a_s = a.as_slice();
    let b_s = b.as_slice();

    for &aa_a in a_s {
        for &aa_b in b_s {
            if aa_a == DEL && aa_b == DEL {
                // both in-frame deletions, no penalty
            } else if aa_a == DEL || aa_b == DEL {
                total_score += DEL_PEN;
            } else if aa_a == FS || aa_b == FS {
                total_score += FS_PEN;
            } else {
                total_score += blosum62_lookup(aa_a, aa_b) as f64;
            }
            total_n += 1;
        }
    }
    if total_n == 0 { 0.0 } else { total_score / total_n as f64 }
}

// ---------------------------------------------------------------------------
// T5: Full-window score computation (IUPAC + BLOSUM62) using flat arrays
// ---------------------------------------------------------------------------

/// Compute the full IUPAC + BLOSUM62 score for a window.
/// `my_notations` and `other_notations` are the flat nucleotide arrays.
/// `other_aas` is precomputed inline AAs for the other side.
/// Returns (total_score, per_codon_iupac, per_codon_blosum) for incremental updates.
pub fn compute_full_score(
    my_notations: &[u8],
    other_notations: &[u8],
    other_aas: &[InlineAA],
    base_score: f64,
) -> (f64, Vec<f64>, Vec<f64>) {
    let n_codons = (my_notations.len() + 2) / 3;
    let mut codon_iupac = Vec::with_capacity(n_codons);
    let mut codon_blosum = Vec::with_capacity(n_codons);
    let mut total = base_score;

    // Process codon by codon
    for ci in 0..n_codons {
        let start = ci * 3;
        let end = (start + 3).min(my_notations.len());

        // IUPAC for this codon's positions
        let mut iupac_sum = 0.0f64;
        for i in start..end {
            if i < other_notations.len() {
                iupac_sum += iupac_score(my_notations[i], other_notations[i]);
            }
        }
        codon_iupac.push(iupac_sum);
        total += iupac_sum;

        // BLOSUM62 for this codon
        let my_aa = translate_codon_inline(&my_notations[start..end]);
        let blosum = if ci < other_aas.len() {
            blosum62_score_inline(&my_aa, &other_aas[ci])
        } else {
            0.0
        };
        codon_blosum.push(blosum);
        total += blosum;
    }

    (total, codon_iupac, codon_blosum)
}

/// Recompute score for a single codon index after gap position change.
/// Returns (new_iupac, new_blosum) for that codon.
#[inline]
pub fn recompute_codon_score(
    my_notations: &[u8],
    other_notations: &[u8],
    other_aa: &InlineAA,
    codon_idx: usize,
) -> (f64, f64) {
    let start = codon_idx * 3;
    let end = (start + 3).min(my_notations.len());

    let mut iupac_sum = 0.0f64;
    for i in start..end {
        if i < other_notations.len() {
            iupac_sum += iupac_score(my_notations[i], other_notations[i]);
        }
    }

    let my_aa = translate_codon_inline(&my_notations[start..end]);
    let blosum = blosum62_score_inline(&my_aa, other_aa);

    (iupac_sum, blosum)
}
