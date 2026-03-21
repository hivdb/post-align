//! Core codon alignment algorithm — Rust port of codon_alignment.py
//! with Option D improvements (D1: precomputed AAs, D4: center-expand).

use std::collections::HashMap;
use std::collections::HashSet;
use std::ops::Range;

use crate::scoring;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A single nucleotide-acid position in an alignment.
///
/// Mirrors `NAPosition` on the Python side.  `notation` stores the
/// ASCII byte of the IUPAC character (or `b'-'`/`b'.'` for gaps),
/// `pos` is the 1-based reference position (`-1` for gaps), and
/// `flag` carries `PositionFlag` metadata.
#[derive(Clone, Copy, Debug)]
pub struct NaPos {
    pub notation: u8,
    pub pos: i32,
    pub flag: i32,
}

impl NaPos {
    #[inline]
    pub fn is_gap(&self) -> bool {
        self.notation == b'-' || self.notation == b'.'
    }
}

/// Per-gap-type scoring adjustments keyed by `(position, gap_length)`.
///
/// Parsed from the Python `gap_placement_score` dictionary in `lib.rs`.
pub struct GapPlacementScore {
    pub refgap: HashMap<(i32, i32), i32>,
    pub seqgap: HashMap<(i32, i32), i32>,
}

const NOGAP: u8 = 0b00;
const REFGAP: u8 = 0b01;
const SEQGAP: u8 = 0b10;

const LEFT: u8 = 0;
const RIGHT: u8 = 1;

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

#[inline]
fn count_gaps(nas: &[NaPos]) -> usize {
    nas.iter().filter(|n| n.is_gap()).count()
}

#[inline]
fn any_has_gap(nas: &[NaPos]) -> bool {
    nas.iter().any(|n| n.is_gap())
}

#[inline]
fn all_have_gap(nas: &[NaPos]) -> bool {
    nas.iter().all(|n| n.is_gap())
}

fn find_first_gap(nas: &[NaPos]) -> i32 {
    for (idx, na) in nas.iter().enumerate() {
        if na.is_gap() {
            return idx as i32;
        }
    }
    -1
}

fn separate_gaps_from_nas(nas: &[NaPos]) -> (Vec<NaPos>, Vec<NaPos>) {
    let mut nongaps = Vec::with_capacity(nas.len());
    let mut gaps = Vec::new();
    for &na in nas {
        if na.is_gap() {
            gaps.push(na);
        } else {
            nongaps.push(na);
        }
    }
    (nongaps, gaps)
}

fn move_gaps_to_center(nas: &[NaPos]) -> Vec<NaPos> {
    let (mut nongaps, gaps) = separate_gaps_from_nas(nas);
    let center = nongaps.len() / 2;
    for (i, g) in gaps.into_iter().enumerate() {
        nongaps.insert(center + i, g);
    }
    nongaps
}

// ---------------------------------------------------------------------------
// Gap window detection
// ---------------------------------------------------------------------------

fn find_windows_with_gap(
    refnas: &[NaPos],
    seqnas: &[NaPos],
    min_gap_distance: i32,
) -> Vec<Range<usize>> {
    let mut first_gap_idx: i64 = -1;
    let mut last_gap_idx: i64 = -1;
    let mut windows: Vec<Range<usize>> = Vec::new();
    for (idx, (refna, seqna)) in refnas.iter().zip(seqnas.iter()).enumerate() {
        if !refna.is_gap() && !seqna.is_gap() {
            continue;
        }
        let idx_i64 = idx as i64;
        if first_gap_idx == -1 {
            first_gap_idx = idx_i64;
            last_gap_idx = idx_i64;
        } else if idx_i64 - last_gap_idx > min_gap_distance as i64 {
            windows.push(first_gap_idx as usize..(last_gap_idx + 1) as usize);
            first_gap_idx = idx_i64;
            last_gap_idx = idx_i64;
        } else {
            last_gap_idx = idx_i64;
        }
    }
    if first_gap_idx > -1 {
        windows.push(first_gap_idx as usize..(last_gap_idx + 1) as usize);
    }
    windows
}

// ---------------------------------------------------------------------------
// Redundant gap removal
// ---------------------------------------------------------------------------

fn remove_n_gaps(nas: &[NaPos], mut n_gaps: usize) -> Vec<NaPos> {
    let mut result = Vec::with_capacity(nas.len());
    for &na in nas {
        if na.is_gap() && n_gaps > 0 {
            n_gaps -= 1;
        } else {
            result.push(na);
        }
    }
    result
}

fn remove_redundant_gaps(
    refnas: &[NaPos],
    seqnas: &[NaPos],
) -> (Vec<NaPos>, Vec<NaPos>) {
    let n = count_gaps(refnas).min(count_gaps(seqnas));
    if n > 0 {
        (remove_n_gaps(refnas, n), remove_n_gaps(seqnas, n))
    } else {
        (refnas.to_vec(), seqnas.to_vec())
    }
}

// ---------------------------------------------------------------------------
// Codon grouping
// ---------------------------------------------------------------------------

/// Split paired NA sequences into codon-aligned groups.
///
/// A new codon group starts every time a non-gap reference position
/// lands on base-pair 0 (mod 3).  Both the reference and sequence
/// sides are split at the same boundaries.
pub fn group_by_codons(
    refnas: &[NaPos],
    seqnas: &[NaPos],
) -> (Vec<Vec<NaPos>>, Vec<Vec<NaPos>>) {
    let mut refcodons: Vec<Vec<NaPos>> = Vec::new();
    let mut seqcodons: Vec<Vec<NaPos>> = Vec::new();
    let mut bp: i32 = -1;
    let mut have_codon = false;
    for (&refna, &seqna) in refnas.iter().zip(seqnas.iter()) {
        if !refna.is_gap() {
            bp = (bp + 1) % 3;
            if bp == 0 {
                refcodons.push(Vec::new());
                seqcodons.push(Vec::new());
                have_codon = true;
            }
        }
        if have_codon {
            refcodons.last_mut().unwrap().push(refna);
            seqcodons.last_mut().unwrap().push(seqna);
        }
    }
    (refcodons, seqcodons)
}

fn find_codon_trim_slice(codons: &[Vec<NaPos>]) -> Range<usize> {
    let mut left_trim: usize = 0;
    for (idx, codon) in codons.iter().enumerate() {
        if all_have_gap(codon) {
            left_trim = idx + 1;
        } else {
            break;
        }
    }
    let codons_len = codons.len();
    let mut right_trim = codons_len;
    for (idx, codon) in codons.iter().rev().enumerate() {
        if all_have_gap(codon) {
            right_trim = codons_len - 1 - idx;
        } else {
            break;
        }
    }
    left_trim..right_trim
}

fn codon_pairs_group_key(refcd: &[NaPos], seqcd: &[NaPos]) -> u8 {
    if any_has_gap(refcd) {
        REFGAP
    } else if any_has_gap(seqcd) {
        SEQGAP
    } else {
        NOGAP
    }
}

fn move_gap_to_codon_end(codons: &[Vec<NaPos>]) -> Vec<Vec<NaPos>> {
    codons
        .iter()
        .map(|codon| {
            let mut nongaps: Vec<NaPos> = codon.iter().filter(|n| !n.is_gap()).copied().collect();
            let gaps: Vec<NaPos> = codon.iter().filter(|n| n.is_gap()).copied().collect();
            nongaps.extend(gaps);
            nongaps
        })
        .collect()
}

fn extend_codons_until_gap(
    ref_codons: &[Vec<NaPos>],
    seq_codons: &[Vec<NaPos>],
    direction: u8,
) -> (Vec<Vec<NaPos>>, Vec<Vec<NaPos>>, usize) {
    let (rcs, scs): (Vec<_>, Vec<_>) = if direction == LEFT {
        (ref_codons.iter().rev().collect(), seq_codons.iter().rev().collect())
    } else {
        (ref_codons.iter().collect(), seq_codons.iter().collect())
    };

    let mut endidx = rcs.len();
    for (idx, (rc, sc)) in rcs.iter().zip(scs.iter()).enumerate() {
        let has_gap = rc.iter().chain(sc.iter()).any(|n| n.is_gap());
        if has_gap {
            endidx = idx;
            break;
        }
    }

    let mut out_ref: Vec<Vec<NaPos>> = rcs[..endidx].iter().map(|c| (*c).clone()).collect();
    let mut out_seq: Vec<Vec<NaPos>> = scs[..endidx].iter().map(|c| (*c).clone()).collect();

    if direction == LEFT {
        out_ref.reverse();
        out_seq.reverse();
    }
    let count = out_ref.len();
    (out_ref, out_seq, count)
}

fn flatten_codons(codons: &[Vec<NaPos>]) -> Vec<NaPos> {
    codons.iter().flat_map(|c| c.iter().copied()).collect()
}

// ---------------------------------------------------------------------------
// D4: Center-expand search order
// ---------------------------------------------------------------------------

fn center_expand_positions(
    orig_gapidx: i32,
    scanstart: usize,
    mynas_len: usize,
    step: usize,
) -> Vec<usize> {
    let mut positions = Vec::new();

    // Snap orig_gapidx to nearest valid step position
    let center = if orig_gapidx >= 0 {
        let raw = orig_gapidx as usize;
        // Round to nearest multiple of step >= scanstart
        let aligned = if raw < scanstart {
            scanstart
        } else {
            ((raw - scanstart) / step) * step + scanstart
        };
        aligned.min(mynas_len)
    } else {
        scanstart
    };

    let mut left = center as i64;
    let mut right = center as i64 + step as i64;
    let mut added_center = false;

    while !added_center || left >= scanstart as i64 || (right as usize) <= mynas_len {
        if !added_center {
            if center <= mynas_len {
                positions.push(center);
            }
            added_center = true;
        }

        if left >= scanstart as i64 && (left as usize) != center {
            if left as usize <= mynas_len {
                positions.push(left as usize);
            }
        }
        if (right as usize) <= mynas_len && (right as usize) != center {
            positions.push(right as usize);
        }

        left -= step as i64;
        right += step as i64;

        if left < scanstart as i64 && right as usize > mynas_len {
            break;
        }
    }
    positions
}

// ---------------------------------------------------------------------------
// Scoring / best-match search (T4 — kept for backward compatibility)
// ---------------------------------------------------------------------------

fn find_best_matches(
    mynas_in: &[NaPos],
    othernas: &[NaPos],
    bp1_indices: &HashSet<usize>,
    gap_type: u8,
    gap_placement_score: &HashMap<(i32, i32), i32>,
    is_start: bool,
    is_end: bool,
) -> Vec<NaPos> {
    let orig_gapidx = find_first_gap(mynas_in);
    let (mynas, mygap) = separate_gaps_from_nas(mynas_in);
    let gaplen = mygap.len();
    let scanstart: usize = if gap_type == REFGAP { 3 } else { 0 };
    let mynas_len = mynas.len();

    // D1: precompute other-side amino acids
    let other_bytes: Vec<u8> = othernas.iter().map(|n| n.notation).collect();
    let other_aas = scoring::translate_codons(&other_bytes);

    // D4: center-expand search order
    let positions = center_expand_positions(orig_gapidx, scanstart, mynas_len, 3);

    let mut max_score: Option<(f64, i32, i32)> = None;
    let mut best_mynas: Option<Vec<NaPos>> = None;

    for idx in positions {
        let mut test_mynas = mynas.clone();
        for (i, &g) in mygap.iter().enumerate() {
            test_mynas.insert(idx + i, g);
        }

        let mut base_score: f64 = -(gaplen as f64);
        if is_start && idx == 0 {
            base_score = 0.0;
        } else if is_end && idx + 3 > mynas_len {
            base_score = 0.0;
        }

        let mut score_val =
            scoring::calc_match_score_precomputed(&test_mynas, othernas, &other_aas, base_score);

        let napos: i32 = if gap_type == REFGAP {
            mynas[idx - 1].pos
        } else {
            othernas[idx].pos
        };

        if let Some(&bonus) = gap_placement_score.get(&(napos, gaplen as i32)) {
            score_val += bonus as f64;
        } else if let Some(&bonus) = gap_placement_score.get(&(napos, 0)) {
            score_val += bonus as f64;
        }

        // Tie-breaking: bp1 bonus (+1, priority 2), orig_gapidx (priority 1),
        // else priority 0. Third element: -idx for leftmost preference.
        let score = if bp1_indices.contains(&idx) {
            (score_val + 1.0, 2i32, -(idx as i32))
        } else if idx as i32 == orig_gapidx {
            (score_val, 1i32, -(idx as i32))
        } else {
            (score_val, 0i32, -(idx as i32))
        };

        let dominated = match &max_score {
            None => false,
            Some(ms) => {
                score.0 < ms.0
                    || (score.0 == ms.0 && score.1 < ms.1)
                    || (score.0 == ms.0 && score.1 == ms.1 && score.2 <= ms.2)
            }
        };
        if !dominated {
            max_score = Some(score);
            best_mynas = Some(test_mynas);
        }
    }

    best_mynas.unwrap_or_else(|| mynas.to_vec())
}

fn paired_find_best_matches(
    refnas: &[NaPos],
    seqnas: &[NaPos],
    gap_type: u8,
    gps: &GapPlacementScore,
    is_seq_start: bool,
    is_seq_end: bool,
) -> (Vec<NaPos>, Vec<NaPos>) {
    // Compute bp1_indices from refnas
    let mut bp1_indices = HashSet::new();
    let mut bp: i32 = 0;
    for (idx, na) in refnas.iter().enumerate() {
        if na.is_gap() {
            continue;
        }
        bp = (bp + 1) % 3;
        if bp == 1 {
            bp1_indices.insert(idx);
        }
    }

    if gap_type == REFGAP {
        let new_refnas =
            find_best_matches(refnas, seqnas, &bp1_indices, REFGAP, &gps.refgap, false, false);
        (new_refnas, seqnas.to_vec())
    } else {
        let new_seqnas = find_best_matches(
            seqnas,
            refnas,
            &bp1_indices,
            SEQGAP,
            &gps.seqgap,
            is_seq_start,
            is_seq_end,
        );
        (refnas.to_vec(), new_seqnas)
    }
}

// ---------------------------------------------------------------------------
// Optimized scoring / best-match search
// ---------------------------------------------------------------------------

/// Build SoA work buffer by inserting gap at position `idx` in nongap arrays.
#[inline]
fn build_work_notation(
    buf: &mut Vec<u8>,
    ng: &[u8],
    gap: &[u8],
    idx: usize,
) {
    buf.clear();
    buf.extend_from_slice(&ng[..idx]);
    buf.extend_from_slice(gap);
    buf.extend_from_slice(&ng[idx..]);
}

/// Optimized find_best_matches: pre-allocated work buffers, InlineAA,
/// precomputed IUPAC table, per-codon incremental score cache.
fn find_best_matches_optimized(
    mynas_in: &[NaPos],
    othernas: &[NaPos],
    bp1_indices: &HashSet<usize>,
    gap_type: u8,
    gap_placement_score: &HashMap<(i32, i32), i32>,
    is_start: bool,
    is_end: bool,
) -> Vec<NaPos> {
    let orig_gapidx = find_first_gap(mynas_in);
    let (mynas_nongap, mygap) = separate_gaps_from_nas(mynas_in);
    let gaplen = mygap.len();
    let scanstart: usize = if gap_type == REFGAP { 3 } else { 0 };
    let mynas_len = mynas_nongap.len();
    let total_len = mynas_len + gaplen;

    // D1: precompute other-side amino acids using InlineAA (no heap per codon)
    let other_n: Vec<u8> = othernas.iter().map(|n| n.notation).collect();
    let other_aas = scoring::translate_codons_inline(&other_n);

    // D4: center-expand search order
    let positions = center_expand_positions(orig_gapidx, scanstart, mynas_len, 3);
    if positions.is_empty() {
        return mynas_in.to_vec();
    }

    // Nongap / gap notation arrays (SoA)
    let ng_n: Vec<u8> = mynas_nongap.iter().map(|n| n.notation).collect();
    let ng_p: Vec<i32> = mynas_nongap.iter().map(|n| n.pos).collect();
    let _ng_f: Vec<i32> = mynas_nongap.iter().map(|n| n.flag).collect();
    let gap_notations: Vec<u8> = mygap.iter().map(|n| n.notation).collect();

    // Pre-allocated work buffer (reused across iterations)
    let mut work_n: Vec<u8> = Vec::with_capacity(total_len);

    // Per-codon score caches
    let n_codons = (total_len + 2) / 3;
    let mut codon_iupac: Vec<f64> = vec![0.0; n_codons];
    let mut codon_blosum: Vec<f64> = vec![0.0; n_codons];

    let mut max_score: Option<(f64, i32, i32)> = None;
    let mut best_idx: usize = positions[0];
    let mut prev_idx: Option<usize> = None;

    // Dummy InlineAA for out-of-bounds codon
    let dummy_aa = scoring::InlineAA { data: [0; 4], len: 0 };

    for &idx in &positions {
        // Build work notation buffer (fast memcpy, no heap alloc)
        build_work_notation(&mut work_n, &ng_n, &gap_notations, idx);

        // Base gap penalty
        let base_score: f64 = if is_start && idx == 0 {
            0.0
        } else if is_end && idx + 3 > mynas_len {
            0.0
        } else {
            -(gaplen as f64)
        };

        let score_val;

        if let Some(old_idx) = prev_idx {
            // Incremental: only recompute codons in the changed range
            let min_changed = old_idx.min(idx);
            let max_changed = (old_idx + gaplen).max(idx + gaplen);
            let first_codon = min_changed / 3;
            let last_codon = ((max_changed + 2) / 3).min(n_codons);

            for ci in first_codon..last_codon {
                let oa = if ci < other_aas.len() { &other_aas[ci] } else { &dummy_aa };
                let (iupac, blosum) = scoring::recompute_codon_score(
                    &work_n, &other_n, oa, ci,
                );
                codon_iupac[ci] = iupac;
                codon_blosum[ci] = blosum;
            }

            // Sum all cached codon scores
            let mut total = base_score;
            for ci in 0..n_codons {
                total += codon_iupac[ci] + codon_blosum[ci];
            }
            score_val = total;
        } else {
            // First position: compute full score and populate caches
            let (total, ci_vec, cb_vec) = scoring::compute_full_score(
                &work_n, &other_n, &other_aas, base_score,
            );
            codon_iupac[..ci_vec.len()].copy_from_slice(&ci_vec);
            codon_blosum[..cb_vec.len()].copy_from_slice(&cb_vec);
            score_val = total;
        }

        prev_idx = Some(idx);

        // Gap placement bonus
        let mut final_score = score_val;
        let napos: i32 = if gap_type == REFGAP {
            if idx > 0 { ng_p[idx - 1] } else { 0 }
        } else {
            if idx < othernas.len() { othernas[idx].pos } else { 0 }
        };

        if let Some(&bonus) = gap_placement_score.get(&(napos, gaplen as i32)) {
            final_score += bonus as f64;
        } else if let Some(&bonus) = gap_placement_score.get(&(napos, 0)) {
            final_score += bonus as f64;
        }

        // Tie-breaking: bp1 bonus (+1, priority 2), orig_gapidx (priority 1)
        let score = if bp1_indices.contains(&idx) {
            (final_score + 1.0, 2i32, -(idx as i32))
        } else if idx as i32 == orig_gapidx {
            (final_score, 1i32, -(idx as i32))
        } else {
            (final_score, 0i32, -(idx as i32))
        };

        let dominated = match &max_score {
            None => false,
            Some(ms) => {
                score.0 < ms.0
                    || (score.0 == ms.0 && score.1 < ms.1)
                    || (score.0 == ms.0 && score.1 == ms.1 && score.2 <= ms.2)
            }
        };
        if !dominated {
            max_score = Some(score);
            best_idx = idx;
        }
    }

    // Reconstruct NaPos result at best_idx
    let mut result = Vec::with_capacity(total_len);
    result.extend_from_slice(&mynas_nongap[..best_idx]);
    result.extend_from_slice(&mygap);
    result.extend_from_slice(&mynas_nongap[best_idx..]);
    result
}

fn paired_find_best_matches_optimized(
    refnas: &[NaPos],
    seqnas: &[NaPos],
    gap_type: u8,
    gps: &GapPlacementScore,
    is_seq_start: bool,
    is_seq_end: bool,
) -> (Vec<NaPos>, Vec<NaPos>) {
    let mut bp1_indices = HashSet::new();
    let mut bp: i32 = 0;
    for (idx, na) in refnas.iter().enumerate() {
        if na.is_gap() {
            continue;
        }
        bp = (bp + 1) % 3;
        if bp == 1 {
            bp1_indices.insert(idx);
        }
    }

    if gap_type == REFGAP {
        let new_refnas = find_best_matches_optimized(
            refnas, seqnas, &bp1_indices, REFGAP, &gps.refgap, false, false,
        );
        (new_refnas, seqnas.to_vec())
    } else {
        let new_seqnas = find_best_matches_optimized(
            seqnas, refnas, &bp1_indices, SEQGAP, &gps.seqgap,
            is_seq_start, is_seq_end,
        );
        (refnas.to_vec(), new_seqnas)
    }
}

// ---------------------------------------------------------------------------
// gather_gaps
// ---------------------------------------------------------------------------

fn gather_gaps(
    mut refnas: Vec<NaPos>,
    mut seqnas: Vec<NaPos>,
    min_gap_distance: i32,
) -> (Vec<NaPos>, Vec<NaPos>) {
    let windows = find_windows_with_gap(&refnas, &seqnas, min_gap_distance);
    // Reverse so replacements don't shift indices
    for window in windows.into_iter().rev() {
        let win_ref = &refnas[window.clone()];
        let win_seq = &seqnas[window.clone()];
        let (wr, ws) = remove_redundant_gaps(win_ref, win_seq);
        let wr = move_gaps_to_center(&wr);
        let ws = move_gaps_to_center(&ws);
        refnas.splice(window.clone(), wr);
        seqnas.splice(window, ws);
    }
    (refnas, seqnas)
}

// ---------------------------------------------------------------------------
// adjust_gap_placement
// ---------------------------------------------------------------------------

fn adjust_gap_placement(
    mut refcodons: Vec<Vec<NaPos>>,
    mut seqcodons: Vec<Vec<NaPos>>,
    window_size: usize,
    gps: &GapPlacementScore,
    is_seq_start: bool,
    is_seq_end: bool,
) -> (Vec<Vec<NaPos>>, Vec<Vec<NaPos>>) {
    let trim_range = find_codon_trim_slice(&seqcodons);

    // Build gap groups: consecutive codons with same gap type
    struct GapGroup {
        gap_type: u8,
        start: usize, // index into refcodons/seqcodons
        end: usize,   // exclusive
    }

    let mut groups: Vec<GapGroup> = Vec::new();

    if trim_range.start < trim_range.end {
        let mut cur_type = NOGAP;
        let mut cur_start = trim_range.start;

        for i in trim_range.clone() {
            let gt = codon_pairs_group_key(&refcodons[i], &seqcodons[i]);
            if i == trim_range.start {
                cur_type = gt;
                cur_start = i;
            } else if gt != cur_type {
                if cur_type != NOGAP {
                    groups.push(GapGroup {
                        gap_type: cur_type,
                        start: cur_start,
                        end: i,
                    });
                }
                cur_type = gt;
                cur_start = i;
            }
        }
        if cur_type != NOGAP {
            groups.push(GapGroup {
                gap_type: cur_type,
                start: cur_start,
                end: trim_range.end,
            });
        }
    }

    // Process each gap group
    for group in &groups {
        let mut start = group.start;
        let mut end = group.end;

        let mut refcds: Vec<Vec<NaPos>> = refcodons[start..end].to_vec();
        let mut seqcds: Vec<Vec<NaPos>> = seqcodons[start..end].to_vec();

        // Extend left
        let left_start = if start > window_size { start - window_size } else { 0 };
        let (ext_r, ext_s, offset) =
            extend_codons_until_gap(&refcodons[left_start..start], &seqcodons[left_start..start], LEFT);
        if offset > 0 {
            let mut new_r = ext_r;
            new_r.extend(refcds);
            refcds = new_r;
            let mut new_s = ext_s;
            new_s.extend(seqcds);
            seqcds = new_s;
            start -= offset;
        }

        // Extend right
        let right_end = (end + window_size).min(refcodons.len());
        let (ext_r, ext_s, offset) =
            extend_codons_until_gap(&refcodons[end..right_end], &seqcodons[end..right_end], RIGHT);
        if offset > 0 {
            refcds.extend(ext_r);
            seqcds.extend(ext_s);
            end += offset;
        }

        let win_refnas = flatten_codons(&refcds);
        let win_seqnas = flatten_codons(&seqcds);

        let (new_refnas, new_seqnas) = paired_find_best_matches(
            &win_refnas,
            &win_seqnas,
            group.gap_type,
            gps,
            is_seq_start && start == trim_range.start,
            is_seq_end && end == trim_range.end,
        );

        let (win_refcodons, win_seqcodons) = group_by_codons(&new_refnas, &new_seqnas);

        // Replace the range in refcodons/seqcodons
        let new_len = win_refcodons.len();
        refcodons.splice(start..end, win_refcodons);
        seqcodons.splice(start..end, win_seqcodons);

        // Note: since we process groups in order and splice may change lengths,
        // subsequent groups' indices might be off. However, the Python original
        // also processes groups in order from the groupby iterator, so the
        // behavior is identical — each group uses the *current* state of
        // refcodons/seqcodons including modifications from previous groups.
        let _ = new_len; // suppress unused warning
    }

    (refcodons, seqcodons)
}

// ---------------------------------------------------------------------------
// realign_gaps — the main entry point called from Python (T4)
// ---------------------------------------------------------------------------

/// Codon-aware gap realignment (baseline Rust port).
///
/// 1. `gather_gaps` — merge nearby gaps and centre them.
/// 2. `group_by_codons` — split into codon-aligned groups.
/// 3. `adjust_gap_placement` — slide each gap group to the
///    position that maximises the IUPAC + BLOSUM62 score.
/// 4. Finalise by moving gaps to codon ends.
pub fn realign_gaps(
    refnas: Vec<NaPos>,
    seqnas: Vec<NaPos>,
    min_gap_distance: i32,
    window_size: i32,
    gps: &GapPlacementScore,
    is_seq_start: bool,
    is_seq_end: bool,
) -> (Vec<NaPos>, Vec<NaPos>) {
    let (refnas, seqnas) = gather_gaps(refnas, seqnas, min_gap_distance);

    let (refcodons, seqcodons) = group_by_codons(&refnas, &seqnas);
    let (refcodons, seqcodons) = adjust_gap_placement(
        refcodons,
        seqcodons,
        window_size as usize,
        gps,
        is_seq_start,
        is_seq_end,
    );

    // Move gaps in seqcodons to codon ends
    let seqcodons = move_gap_to_codon_end(&seqcodons);

    (flatten_codons(&refcodons), flatten_codons(&seqcodons))
}

// ---------------------------------------------------------------------------
// adjust_gap_placement using optimized scoring
// ---------------------------------------------------------------------------

fn adjust_gap_placement_optimized(
    mut refcodons: Vec<Vec<NaPos>>,
    mut seqcodons: Vec<Vec<NaPos>>,
    window_size: usize,
    gps: &GapPlacementScore,
    is_seq_start: bool,
    is_seq_end: bool,
) -> (Vec<Vec<NaPos>>, Vec<Vec<NaPos>>) {
    let trim_range = find_codon_trim_slice(&seqcodons);

    struct GapGroup {
        gap_type: u8,
        start: usize,
        end: usize,
    }

    let mut groups: Vec<GapGroup> = Vec::new();

    if trim_range.start < trim_range.end {
        let mut cur_type = NOGAP;
        let mut cur_start = trim_range.start;

        for i in trim_range.clone() {
            let gt = codon_pairs_group_key(&refcodons[i], &seqcodons[i]);
            if i == trim_range.start {
                cur_type = gt;
                cur_start = i;
            } else if gt != cur_type {
                if cur_type != NOGAP {
                    groups.push(GapGroup {
                        gap_type: cur_type,
                        start: cur_start,
                        end: i,
                    });
                }
                cur_type = gt;
                cur_start = i;
            }
        }
        if cur_type != NOGAP {
            groups.push(GapGroup {
                gap_type: cur_type,
                start: cur_start,
                end: trim_range.end,
            });
        }
    }

    for group in &groups {
        let mut start = group.start;
        let mut end = group.end;

        let mut refcds: Vec<Vec<NaPos>> = refcodons[start..end].to_vec();
        let mut seqcds: Vec<Vec<NaPos>> = seqcodons[start..end].to_vec();

        // Extend left
        let left_start = if start > window_size { start - window_size } else { 0 };
        let (ext_r, ext_s, offset) =
            extend_codons_until_gap(&refcodons[left_start..start], &seqcodons[left_start..start], LEFT);
        if offset > 0 {
            let mut new_r = ext_r;
            new_r.extend(refcds);
            refcds = new_r;
            let mut new_s = ext_s;
            new_s.extend(seqcds);
            seqcds = new_s;
            start -= offset;
        }

        // Extend right
        let right_end = (end + window_size).min(refcodons.len());
        let (ext_r, ext_s, offset) =
            extend_codons_until_gap(&refcodons[end..right_end], &seqcodons[end..right_end], RIGHT);
        if offset > 0 {
            refcds.extend(ext_r);
            seqcds.extend(ext_s);
            end += offset;
        }

        let win_refnas = flatten_codons(&refcds);
        let win_seqnas = flatten_codons(&seqcds);

        let (new_refnas, new_seqnas) = paired_find_best_matches_optimized(
            &win_refnas,
            &win_seqnas,
            group.gap_type,
            gps,
            is_seq_start && start == trim_range.start,
            is_seq_end && end == trim_range.end,
        );

        let (win_refcodons, win_seqcodons) = group_by_codons(&new_refnas, &new_seqnas);

        let new_len = win_refcodons.len();
        refcodons.splice(start..end, win_refcodons);
        seqcodons.splice(start..end, win_seqcodons);
        let _ = new_len;
    }

    (refcodons, seqcodons)
}

// ---------------------------------------------------------------------------
// realign_gaps_optimized — optimized entry point
// ---------------------------------------------------------------------------

/// Optimised codon-aware gap realignment.
///
/// Same algorithm as [`realign_gaps`] but uses the optimised scoring
/// path: precomputed IUPAC lookup table, `InlineAA` codon translation
/// (no heap allocation per codon), incremental per-codon score cache,
/// and centre-expand search order.
pub fn realign_gaps_optimized(
    refnas: Vec<NaPos>,
    seqnas: Vec<NaPos>,
    min_gap_distance: i32,
    window_size: i32,
    gps: &GapPlacementScore,
    is_seq_start: bool,
    is_seq_end: bool,
) -> (Vec<NaPos>, Vec<NaPos>) {
    let (refnas, seqnas) = gather_gaps(refnas, seqnas, min_gap_distance);

    let (refcodons, seqcodons) = group_by_codons(&refnas, &seqnas);
    let (refcodons, seqcodons) = adjust_gap_placement_optimized(
        refcodons,
        seqcodons,
        window_size as usize,
        gps,
        is_seq_start,
        is_seq_end,
    );

    let seqcodons = move_gap_to_codon_end(&seqcodons);

    (flatten_codons(&refcodons), flatten_codons(&seqcodons))
}

