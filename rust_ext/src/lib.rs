//! postalign_rs — Rust core for codon alignment via PyO3.

mod align;
mod scoring;

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use std::collections::HashMap;

use align::{GapPlacementScore, NaPos};
use rayon::prelude::*;

/// Convert Python gap_placement_score dict to Rust struct.
///
/// Python shape: {gap_type_int: {(napos, gaplen): score, ...}, ...}
/// where gap_type 1 = REFGAP, 2 = SEQGAP.
fn parse_gap_placement_score(py_dict: &Bound<'_, PyDict>) -> PyResult<GapPlacementScore> {
    let mut refgap = HashMap::new();
    let mut seqgap = HashMap::new();

    for (key, value) in py_dict.iter() {
        let gap_type: i32 = key.extract()?;
        let inner: &Bound<'_, PyDict> = value.downcast()?;
        let target = if gap_type == 1 { &mut refgap } else { &mut seqgap };
        for (k2, v2) in inner.iter() {
            let tup: &Bound<'_, PyTuple> = k2.downcast()?;
            let napos: i32 = tup.get_item(0)?.extract()?;
            let gaplen: i32 = tup.get_item(1)?.extract()?;
            let score: i32 = v2.extract()?;
            target.insert((napos, gaplen), score);
        }
    }

    Ok(GapPlacementScore { refgap, seqgap })
}

/// Codon-aware gap realignment (baseline Rust port).
///
/// Accepts flat arrays decomposed from `NAPosition` lists, runs
/// `gather_gaps` → `group_by_codons` → `adjust_gap_placement`, and
/// returns the realigned flat arrays as a 6-element Python tuple:
/// `(ref_notations, ref_positions, ref_flags,
///   seq_notations, seq_positions, seq_flags)`.
#[pyfunction]
fn realign_gaps(
    py: Python<'_>,
    ref_notations: Vec<u8>,
    ref_positions: Vec<i32>,
    ref_flags: Vec<i32>,
    seq_notations: Vec<u8>,
    seq_positions: Vec<i32>,
    seq_flags: Vec<i32>,
    min_gap_distance: i32,
    window_size: i32,
    gap_placement_score: &Bound<'_, PyDict>,
    is_seq_start: bool,
    is_seq_end: bool,
) -> PyResult<PyObject> {
    let n_ref = ref_notations.len();
    let n_seq = seq_notations.len();

    // Build NaPos vectors
    let mut refnas = Vec::with_capacity(n_ref);
    for i in 0..n_ref {
        refnas.push(NaPos {
            notation: ref_notations[i],
            pos: ref_positions[i],
            flag: ref_flags[i],
        });
    }
    let mut seqnas = Vec::with_capacity(n_seq);
    for i in 0..n_seq {
        seqnas.push(NaPos {
            notation: seq_notations[i],
            pos: seq_positions[i],
            flag: seq_flags[i],
        });
    }

    let gps = parse_gap_placement_score(gap_placement_score)?;

    let (out_ref, out_seq) = align::realign_gaps(
        refnas,
        seqnas,
        min_gap_distance,
        window_size,
        &gps,
        is_seq_start,
        is_seq_end,
    );

    // Convert back to Python lists
    let ref_n = PyList::new(py, out_ref.iter().map(|n| n.notation))?;
    let ref_p = PyList::new(py, out_ref.iter().map(|n| n.pos))?;
    let ref_f = PyList::new(py, out_ref.iter().map(|n| n.flag))?;
    let seq_n = PyList::new(py, out_seq.iter().map(|n| n.notation))?;
    let seq_p = PyList::new(py, out_seq.iter().map(|n| n.pos))?;
    let seq_f = PyList::new(py, out_seq.iter().map(|n| n.flag))?;

    Ok(PyTuple::new(py, &[
        ref_n.into_any(),
        ref_p.into_any(),
        ref_f.into_any(),
        seq_n.into_any(),
        seq_p.into_any(),
        seq_f.into_any(),
    ])?.into())
}

/// Optimised codon-aware gap realignment.
///
/// Drop-in replacement for [`realign_gaps`] that uses the optimised
/// scoring path: compile-time IUPAC lookup table, `InlineAA`
/// stack-allocated codon translation, incremental per-codon score
/// cache, and centre-expand search order.
#[pyfunction]
fn realign_gaps_optimized(
    py: Python<'_>,
    ref_notations: Vec<u8>,
    ref_positions: Vec<i32>,
    ref_flags: Vec<i32>,
    seq_notations: Vec<u8>,
    seq_positions: Vec<i32>,
    seq_flags: Vec<i32>,
    min_gap_distance: i32,
    window_size: i32,
    gap_placement_score: &Bound<'_, PyDict>,
    is_seq_start: bool,
    is_seq_end: bool,
) -> PyResult<PyObject> {
    let n_ref = ref_notations.len();
    let n_seq = seq_notations.len();

    let mut refnas = Vec::with_capacity(n_ref);
    for i in 0..n_ref {
        refnas.push(NaPos {
            notation: ref_notations[i],
            pos: ref_positions[i],
            flag: ref_flags[i],
        });
    }
    let mut seqnas = Vec::with_capacity(n_seq);
    for i in 0..n_seq {
        seqnas.push(NaPos {
            notation: seq_notations[i],
            pos: seq_positions[i],
            flag: seq_flags[i],
        });
    }

    let gps = parse_gap_placement_score(gap_placement_score)?;

    let (out_ref, out_seq) = align::realign_gaps_optimized(
        refnas,
        seqnas,
        min_gap_distance,
        window_size,
        &gps,
        is_seq_start,
        is_seq_end,
    );

    let ref_n = PyList::new(py, out_ref.iter().map(|n| n.notation))?;
    let ref_p = PyList::new(py, out_ref.iter().map(|n| n.pos))?;
    let ref_f = PyList::new(py, out_ref.iter().map(|n| n.flag))?;
    let seq_n = PyList::new(py, out_seq.iter().map(|n| n.notation))?;
    let seq_p = PyList::new(py, out_seq.iter().map(|n| n.pos))?;
    let seq_f = PyList::new(py, out_seq.iter().map(|n| n.flag))?;

    Ok(PyTuple::new(py, &[
        ref_n.into_any(),
        ref_p.into_any(),
        ref_f.into_any(),
        seq_n.into_any(),
        seq_p.into_any(),
        seq_f.into_any(),
    ])?.into())
}

// ---------------------------------------------------------------------------
// Phase C: Full codon_align in Rust (boundary detection + realign_gaps)
// ---------------------------------------------------------------------------

/// Helper: build NaPos vec from flat arrays.
#[inline]
fn build_napos_vec(notations: &[u8], positions: &[i32], flags: &[i32]) -> Vec<NaPos> {
    let n = notations.len();
    let mut v = Vec::with_capacity(n);
    for i in 0..n {
        v.push(NaPos {
            notation: notations[i],
            pos: positions[i],
            flag: flags[i],
        });
    }
    v
}

/// Helper: convert NaPos vec back to Python tuple of 3 lists.
fn napos_to_py(py: Python<'_>, nas: &[NaPos]) -> PyResult<(Py<PyList>, Py<PyList>, Py<PyList>)> {
    let n_list = PyList::new(py, nas.iter().map(|n| n.notation))?;
    let p_list = PyList::new(py, nas.iter().map(|n| n.pos))?;
    let f_list = PyList::new(py, nas.iter().map(|n| n.flag))?;
    Ok((n_list.into(), p_list.into(), f_list.into()))
}

/// Rust port of NAPosition.min_pos
#[inline]
fn min_pos(positions: &[i32]) -> i32 {
    for &p in positions {
        if p > 0 { return p; }
    }
    -1
}

/// Rust port of NAPosition.max_pos
#[inline]
fn max_pos(positions: &[i32]) -> i32 {
    for &p in positions.iter().rev() {
        if p > 0 { return p; }
    }
    -1
}

/// Rust port of NAPosition.min_nongap_index(nas, start)
#[inline]
fn min_nongap_index(positions: &[i32], start: usize) -> i64 {
    for idx in start..positions.len() {
        if positions[idx] > 0 { return idx as i64; }
    }
    -1
}

/// Rust port of NAPosition.max_nongap_index(nas)
#[inline]
fn max_nongap_index(positions: &[i32]) -> i64 {
    for idx in (0..positions.len()).rev() {
        if positions[idx] > 0 { return idx as i64; }
    }
    -1
}

/// Rust port of _pos2index (find first/last index with given pos)
#[inline]
fn pos2index_first(positions: &[i32], pos: i32) -> i64 {
    for (idx, &p) in positions.iter().enumerate() {
        if p == pos { return idx as i64; }
    }
    -1
}

#[inline]
fn pos2index_last(positions: &[i32], pos: i32) -> i64 {
    for idx in (0..positions.len()).rev() {
        if positions[idx] == pos { return idx as i64; }
    }
    -1
}

/// Rust port of _posrange2indexrange with include_boundary_gaps=True
fn posrange2indexrange(
    notations: &[u8], positions: &[i32],
    pos_start: i32, pos_end: i32,
) -> (usize, usize) {
    let mp = max_pos(positions);
    let mnp = min_pos(positions);
    if mp < 0 || mnp < 0 {
        return (0, 0);
    }

    let (mut idx_start, mut idx_end): (i64, i64);

    if pos_start > mp {
        let v = pos2index_last(positions, mp);
        idx_start = v;
        idx_end = v;
    } else if pos_end < mnp {
        let v = pos2index_first(positions, mnp);
        idx_start = v;
        idx_end = v;
    } else {
        let ps = if mnp > pos_start { mnp } else { pos_start };
        let pe = if mp < pos_end { mp } else { pos_end };

        idx_start = -1;
        for pos in ps..=pe {
            idx_start = pos2index_first(positions, pos);
            if idx_start > -1 { break; }
        }

        idx_end = -1;
        for pos in (ps..=pe).rev() {
            let v = pos2index_last(positions, pos);
            if v > -1 { idx_end = v + 1; break; }
        }
    }

    // include_boundary_gaps: extend to cover adjacent gaps
    let mut is = idx_start as usize;
    let mut ie = idx_end as usize;
    let n = notations.len();

    while is > 0 && (notations[is - 1] == b'-' || notations[is - 1] == b'.') {
        is -= 1;
    }
    while ie < n && (notations[ie] == b'-' || notations[ie] == b'.') {
        ie += 1;
    }

    (is, ie)
}

/// Rust port of any_has_gap
#[inline]
fn any_has_gap(notations: &[u8]) -> bool {
    notations.iter().any(|&b| b == b'-' || b == b'.')
}

// ---------------------------------------------------------------------------
// Pure-Rust core: boundary detection + realign (no PyO3 types)
// ---------------------------------------------------------------------------

/// Result of a single codon_align operation.
struct AlignResult {
    idx_start: usize,
    idx_end: usize,
    out_ref: Vec<NaPos>,
    out_seq: Vec<NaPos>,
}

/// Pure-Rust codon_align: boundary detection + optimised realign_gaps.
/// Returns None if no alignment needed.
fn codon_align_core(
    ref_notations: &[u8],
    ref_positions: &[i32],
    ref_flags: &[i32],
    seq_notations: &[u8],
    seq_positions: &[i32],
    seq_flags: &[i32],
    min_gap_distance: i32,
    window_size: i32,
    gps: &GapPlacementScore,
    ref_start: i32,
    ref_end: i32,
) -> Option<AlignResult> {
    let n_seq = seq_notations.len();

    let (ref_idx_start, ref_idx_end) = posrange2indexrange(
        ref_notations, ref_positions, ref_start, ref_end,
    );

    let seq_idx_start: usize = 0;
    let seq_idx_end: usize = n_seq;

    let mut idx_start = if ref_idx_start > seq_idx_start {
        ref_idx_start
    } else {
        seq_idx_start
    };

    loop {
        let test = min_nongap_index(ref_positions, idx_start);
        if test < 0 { break; }
        let ti = test as usize;
        if (ref_positions[ti] - ref_start) % 3 == 0 { break; }
        idx_start = ti + 1;
    }

    let idx_end = if ref_idx_end < seq_idx_end {
        ref_idx_end
    } else {
        seq_idx_end
    };

    if idx_start == idx_end {
        return None;
    }

    let is_seq_start = (idx_start as i64) <= min_nongap_index(seq_positions, 0);
    let is_seq_end = (idx_end as i64) > max_nongap_index(seq_positions);

    let ref_window = &ref_notations[idx_start..idx_end];
    let seq_window = &seq_notations[idx_start..idx_end];
    if !any_has_gap(ref_window) && !any_has_gap(seq_window) {
        return None;
    }

    let refnas = build_napos_vec(
        &ref_notations[idx_start..idx_end],
        &ref_positions[idx_start..idx_end],
        &ref_flags[idx_start..idx_end],
    );
    let seqnas = build_napos_vec(
        &seq_notations[idx_start..idx_end],
        &seq_positions[idx_start..idx_end],
        &seq_flags[idx_start..idx_end],
    );

    let (out_ref, out_seq) = align::realign_gaps_optimized(
        refnas, seqnas,
        min_gap_distance, window_size,
        gps, is_seq_start, is_seq_end,
    );

    Some(AlignResult { idx_start, idx_end, out_ref, out_seq })
}

/// Convert AlignResult to Python tuple, or return None.
fn align_result_to_py(py: Python<'_>, res: Option<AlignResult>) -> PyResult<PyObject> {
    match res {
        None => Ok(py.None()),
        Some(r) => {
            let (rn, rp, rf) = napos_to_py(py, &r.out_ref)?;
            let (sn, sp, sf) = napos_to_py(py, &r.out_seq)?;
            Ok(PyTuple::new(py, &[
                r.idx_start.into_pyobject(py)?.into_any(),
                r.idx_end.into_pyobject(py)?.into_any(),
                rn.bind(py).clone().into_any(),
                rp.bind(py).clone().into_any(),
                rf.bind(py).clone().into_any(),
                sn.bind(py).clone().into_any(),
                sp.bind(py).clone().into_any(),
                sf.bind(py).clone().into_any(),
            ])?.into())
        }
    }
}

/// Full codon alignment in Rust for a single sequence pair.
///
/// Performs boundary detection (equivalent to Python's
/// `_posrange2indexrange`) followed by optimised gap realignment.
/// Returns `None` if no alignment is needed, otherwise an 8-tuple:
/// `(idx_start, idx_end, ref_n, ref_p, ref_f, seq_n, seq_p, seq_f)`.
#[pyfunction]
fn codon_align_full(
    py: Python<'_>,
    ref_notations: Vec<u8>,
    ref_positions: Vec<i32>,
    ref_flags: Vec<i32>,
    seq_notations: Vec<u8>,
    seq_positions: Vec<i32>,
    seq_flags: Vec<i32>,
    min_gap_distance: i32,
    window_size: i32,
    gap_placement_score: &Bound<'_, PyDict>,
    ref_start: i32,
    ref_end: i32,
) -> PyResult<PyObject> {
    let gps = parse_gap_placement_score(gap_placement_score)?;
    let result = codon_align_core(
        &ref_notations, &ref_positions, &ref_flags,
        &seq_notations, &seq_positions, &seq_flags,
        min_gap_distance, window_size, &gps,
        ref_start, ref_end,
    );
    align_result_to_py(py, result)
}

// ---------------------------------------------------------------------------
// Batch API with rayon parallelism
// ---------------------------------------------------------------------------

/// Input data for one sequence pair (owned, Send+Sync safe).
struct BatchInput {
    ref_notations: Vec<u8>,
    ref_positions: Vec<i32>,
    ref_flags: Vec<i32>,
    seq_notations: Vec<u8>,
    seq_positions: Vec<i32>,
    seq_flags: Vec<i32>,
    ref_start: i32,
    ref_end: i32,
}

/// Batch codon alignment with rayon parallelism.
///
/// Processes *N* sequence pairs in parallel, releasing the GIL
/// during the compute phase.  Each input item is an 8-tuple of
/// flat arrays `(ref_n, ref_p, ref_f, seq_n, seq_p, seq_f,
/// ref_start, ref_end)`.  Shared parameters (`min_gap_distance`,
/// `window_size`, `gap_placement_score`) are parsed once.
///
/// Returns a Python list of results, each `None` or an 8-tuple.
#[pyfunction]
fn codon_align_batch(
    py: Python<'_>,
    items: &Bound<'_, PyList>,
    min_gap_distance: i32,
    window_size: i32,
    gap_placement_score: &Bound<'_, PyDict>,
) -> PyResult<PyObject> {
    // 1) Parse shared GPS once
    let gps = parse_gap_placement_score(gap_placement_score)?;

    // 2) Extract all inputs from Python (serial, unavoidable GIL work)
    let n = items.len();
    let mut inputs: Vec<BatchInput> = Vec::with_capacity(n);

    for i in 0..n {
        let item = items.get_item(i)?;
        let tup = item.downcast::<PyTuple>()?;
        inputs.push(BatchInput {
            ref_notations: tup.get_item(0)?.extract()?,
            ref_positions: tup.get_item(1)?.extract()?,
            ref_flags: tup.get_item(2)?.extract()?,
            seq_notations: tup.get_item(3)?.extract()?,
            seq_positions: tup.get_item(4)?.extract()?,
            seq_flags: tup.get_item(5)?.extract()?,
            ref_start: tup.get_item(6)?.extract()?,
            ref_end: tup.get_item(7)?.extract()?,
        });
    }

    // 3) Release GIL and process in parallel
    let results: Vec<Option<AlignResult>> = py.allow_threads(|| {
        inputs.par_iter().map(|inp| {
            codon_align_core(
                &inp.ref_notations, &inp.ref_positions, &inp.ref_flags,
                &inp.seq_notations, &inp.seq_positions, &inp.seq_flags,
                min_gap_distance, window_size, &gps,
                inp.ref_start, inp.ref_end,
            )
        }).collect()
    });

    // 4) Convert results back to Python (serial, needs GIL)
    let out = PyList::empty(py);
    for res in results {
        let py_obj = align_result_to_py(py, res)?;
        out.append(py_obj)?;
    }

    Ok(out.into_any().unbind())
}

/// Python module definition.
#[pymodule]
fn postalign_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(realign_gaps, m)?)?;
    m.add_function(wrap_pyfunction!(realign_gaps_optimized, m)?)?;
    m.add_function(wrap_pyfunction!(codon_align_full, m)?)?;
    m.add_function(wrap_pyfunction!(codon_align_batch, m)?)?;
    Ok(())
}
