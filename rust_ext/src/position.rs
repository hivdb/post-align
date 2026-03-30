//! Rust-backed NAPosition and AAPosition types exposed to Python via PyO3.

use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyList, PySlice};

/// Gap character constants (ASCII dash and dot).
const GAP_DASH: u8 = b'-';
const GAP_DOT: u8 = b'.';

#[inline]
fn is_gap_byte(b: u8) -> bool {
    b == GAP_DASH || b == GAP_DOT
}

/// Clone helper for `Option<PyObject>`.  `Py::clone_ref` requires
/// the GIL token, so we grab it ourselves.
#[inline]
fn clone_opt_pyobj(opt: &Option<PyObject>) -> Option<PyObject> {
    opt.as_ref().map(|o| Python::with_gil(|py| o.clone_ref(py)))
}

// ---------------------------------------------------------------------------
// NAPosition
// ---------------------------------------------------------------------------

/// A single nucleotide-acid position.
///
/// Drop-in replacement for the Cython `NAPosition` class.  All fields,
/// constructor parameters, dunder methods, and static/class methods
/// are preserved so that existing Python code works unchanged.
#[pyclass(module = "postalign_rs")]
#[derive(Debug)]
pub struct NAPosition {
    #[pyo3(get, set)]
    pub notation: i32,
    #[pyo3(get, set)]
    pub pos: i32,
    // NOTE: we do NOT use `#[pyo3(set)]` on `flag` because PyO3
    // would generate `set_flag` which clashes with our static method
    // of the same name.  We implement the setter manually instead.
    #[pyo3(get)]
    pub flag: i32,
    #[pyo3(get)]
    pub is_gap: bool,
    #[pyo3(get, set)]
    pub payload: Option<PyObject>,
}

impl Clone for NAPosition {
    fn clone(&self) -> Self {
        NAPosition {
            notation: self.notation,
            pos: self.pos,
            flag: self.flag,
            is_gap: self.is_gap,
            payload: clone_opt_pyobj(&self.payload),
        }
    }
}

#[pymethods]
impl NAPosition {
    #[new]
    #[pyo3(signature = (notation, pos, flag, payload=None))]
    fn new(notation: i32, pos: i32, flag: i32, payload: Option<PyObject>) -> Self {
        NAPosition {
            notation,
            pos,
            flag,
            is_gap: is_gap_byte(notation as u8),
            payload,
        }
    }

    /// Manual setter for `flag` so that `obj.flag = x` works from
    /// Python without colliding with the `set_flag` static method.
    #[setter(flag)]
    fn _set_flag(&mut self, value: i32) {
        self.flag = value;
    }

    fn __str__(&self) -> String {
        String::from(self.notation as u8 as char)
    }

    fn __bytes__<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &[self.notation as u8])
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        match &self.payload {
            None => Ok(format!(
                "NAPosition({:?}, {:?}, {:?}, None)",
                self.notation, self.pos, self.flag
            )),
            Some(obj) => {
                let repr_str = obj.bind(py).repr()?.to_string();
                Ok(format!(
                    "NAPosition({:?}, {:?}, {:?}, {})",
                    self.notation, self.pos, self.flag, repr_str
                ))
            }
        }
    }

    fn __copy__(&self) -> Self {
        self.clone()
    }

    // -- Class methods -------------------------------------------------------

    /// Create a list of gap NAPositions.
    #[classmethod]
    #[pyo3(signature = (gaplen,))]
    fn init_gaps(_cls: &Bound<'_, pyo3::types::PyType>, gaplen: usize) -> Vec<NAPosition> {
        (0..gaplen)
            .map(|_| NAPosition {
                notation: GAP_DASH as i32,
                pos: -1,
                flag: 0, // PositionFlag.NONE
                is_gap: true,
                payload: None,
            })
            .collect()
    }

    /// Create NAPosition list from a byte sequence.
    #[classmethod]
    #[pyo3(signature = (seq_text, seq_payload=None))]
    fn init_from_bytes(
        _cls: &Bound<'_, pyo3::types::PyType>,
        seq_text: &[u8],
        seq_payload: Option<Vec<PyObject>>,
    ) -> Vec<NAPosition> {
        let upper: Vec<u8> = seq_text.iter().map(|b| b.to_ascii_uppercase()).collect();
        let payloads = seq_payload.unwrap_or_default();
        let mut offset: i32 = 1;
        upper
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                let is_g = is_gap_byte(b);
                let pos = if is_g {
                    -1
                } else {
                    let p = offset;
                    offset += 1;
                    p
                };
                NAPosition {
                    notation: b as i32,
                    pos,
                    flag: 0,
                    is_gap: is_g,
                    payload: payloads.get(i).map(|o| {
                        Python::with_gil(|py| o.clone_ref(py))
                    }),
                }
            })
            .collect()
    }

    // -- Static methods operating on lists -----------------------------------
    //
    // All methods accept `&Bound<'_, PyList>` and borrow each element
    // lazily via `downcast::<NAPosition>()?.borrow()`.  This avoids the
    // O(n) upfront extraction cost of `Vec<PyRef<'_, NAPosition>>`
    // and enables truly early-return behaviour for methods like
    // `min_pos` / `max_pos`.

    /// Return the first positive position, or -1.
    #[staticmethod]
    fn min_pos(nas: &Bound<'_, PyList>) -> PyResult<i32> {
        for i in 0..nas.len() {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            let r = na.borrow();
            if r.pos > 0 {
                return Ok(r.pos);
            }
        }
        Ok(-1)
    }

    /// Return the last positive position, or -1.
    #[staticmethod]
    fn max_pos(nas: &Bound<'_, PyList>) -> PyResult<i32> {
        let len = nas.len();
        for i in (0..len).rev() {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            let r = na.borrow();
            if r.pos > 0 {
                return Ok(r.pos);
            }
        }
        Ok(-1)
    }

    /// Index of first non-gap element in [start..stop), or -1.
    #[staticmethod]
    #[pyo3(signature = (nas, start=-1, stop=-1))]
    fn min_nongap_index(
        nas: &Bound<'_, PyList>,
        start: i32,
        stop: i32,
    ) -> PyResult<i32> {
        let len = nas.len() as i32;
        for idx in 0..len {
            if start > -1 && idx < start {
                continue;
            }
            if stop > -1 && idx >= stop {
                break;
            }
            let na = nas.get_item(idx as usize)?.downcast_into::<NAPosition>()?;
            if na.borrow().pos > 0 {
                return Ok(idx);
            }
        }
        Ok(-1)
    }

    /// Index of last non-gap element in [start..stop), or -1.
    #[staticmethod]
    #[pyo3(signature = (nas, start=-1, stop=-1))]
    fn max_nongap_index(
        nas: &Bound<'_, PyList>,
        start: i32,
        stop: i32,
    ) -> PyResult<i32> {
        let length = nas.len() as i32;
        for revidx in 0..length {
            let idx = length - 1 - revidx;
            if start > -1 && idx < start {
                break;
            }
            if stop > -1 && idx >= stop {
                continue;
            }
            let na = nas.get_item(idx as usize)?.downcast_into::<NAPosition>()?;
            if na.borrow().pos > 0 {
                return Ok(idx);
            }
        }
        Ok(-1)
    }

    /// Set flag (bitwise OR) on every element.
    #[staticmethod]
    #[pyo3(name = "set_flag")]
    fn set_flag_static(nas: &Bound<'_, PyList>, flag: i32) -> PyResult<()> {
        for i in 0..nas.len() {
            let item = nas.get_item(i)?;
            let mut na = item.downcast::<NAPosition>()?.borrow_mut();
            na.flag |= flag;
        }
        Ok(())
    }

    /// True if any element has the given flag bits set.
    #[staticmethod]
    fn any_has_flag(nas: &Bound<'_, PyList>, flag: i32) -> PyResult<bool> {
        for i in 0..nas.len() {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            if (na.borrow().flag & flag) != 0 {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// True if all elements have the given flag bits set.
    #[staticmethod]
    fn all_have_flag(nas: &Bound<'_, PyList>, flag: i32) -> PyResult<bool> {
        for i in 0..nas.len() {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            if (na.borrow().flag & flag) == 0 {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Number of gap elements.
    #[staticmethod]
    fn count_gaps(nas: &Bound<'_, PyList>) -> PyResult<usize> {
        let mut count = 0;
        for i in 0..nas.len() {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            if na.borrow().is_gap {
                count += 1;
            }
        }
        Ok(count)
    }

    /// Number of non-gap elements.
    #[staticmethod]
    fn count_nongaps(nas: &Bound<'_, PyList>) -> PyResult<usize> {
        let mut count = 0;
        for i in 0..nas.len() {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            if !na.borrow().is_gap {
                count += 1;
            }
        }
        Ok(count)
    }

    /// True if any element is a gap.
    #[staticmethod]
    fn any_has_gap(nas: &Bound<'_, PyList>) -> PyResult<bool> {
        for i in 0..nas.len() {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            if na.borrow().is_gap {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// True if all elements are gaps.
    #[staticmethod]
    fn all_have_gap(nas: &Bound<'_, PyList>) -> PyResult<bool> {
        for i in 0..nas.len() {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            if !na.borrow().is_gap {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Return only non-gap elements (as references to existing objects).
    ///
    /// Uses `PyList::get_item_unchecked` + direct downcast pointer
    /// access for maximum throughput on the filter loop.
    #[staticmethod]
    fn remove_gaps<'py>(nas: &Bound<'py, PyList>) -> PyResult<Bound<'py, PyList>> {
        let py = nas.py();
        let len = nas.len();
        let result = PyList::empty(py);
        for i in 0..len {
            // SAFETY: i is in 0..len
            let item = unsafe { nas.get_item_unchecked(i) };
            let na: &Bound<'_, NAPosition> = item.downcast()?;
            if !na.borrow().is_gap {
                result.append(&item)?;
            }
        }
        Ok(result)
    }

    /// Concatenate notations into bytes.
    #[staticmethod]
    fn as_bytes<'py>(
        py: Python<'py>,
        nas: &Bound<'py, PyList>,
    ) -> PyResult<Bound<'py, PyBytes>> {
        let mut v: Vec<u8> = Vec::with_capacity(nas.len());
        for i in 0..nas.len() {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            v.push(na.borrow().notation as u8);
        }
        Ok(PyBytes::new(py, &v))
    }

    /// Concatenate notations into a string.
    #[classmethod]
    fn as_str(
        _cls: &Bound<'_, pyo3::types::PyType>,
        nas: &Bound<'_, PyList>,
    ) -> PyResult<String> {
        let mut s = String::with_capacity(nas.len());
        for i in 0..nas.len() {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            s.push(na.borrow().notation as u8 as char);
        }
        Ok(s)
    }

    /// Map position range to index range.
    #[staticmethod]
    #[pyo3(signature = (nas, pos_start, pos_end, include_boundary_gaps=false))]
    fn posrange2indexrange(
        nas: &Bound<'_, PyList>,
        pos_start: i32,
        pos_end: i32,
        include_boundary_gaps: bool,
    ) -> PyResult<(usize, usize)> {
        let len = nas.len();
        let max_p = Self::_max_pos_list(nas)?;
        let min_p = Self::_min_pos_list(nas)?;
        if max_p < 0 || min_p < 0 {
            return Ok((0, 0));
        }

        let mut idx_start: i32;
        let mut idx_end: i32;

        if pos_start > max_p {
            let v = Self::_pos2index_list(nas, max_p, false)?;
            idx_start = v;
            idx_end = v;
        } else if pos_end < min_p {
            let v = Self::_pos2index_list(nas, min_p, true)?;
            idx_start = v;
            idx_end = v;
        } else {
            let ps = if min_p > pos_start { min_p } else { pos_start };
            let pe = if max_p < pos_end { max_p } else { pos_end };

            idx_start = -1;
            for pos in ps..=pe {
                idx_start = Self::_pos2index_list(nas, pos, true)?;
                if idx_start > -1 {
                    break;
                }
            }

            idx_end = -1;
            for pos in (ps..=pe).rev() {
                let v = Self::_pos2index_list(nas, pos, false)?;
                if v > -1 {
                    idx_end = v + 1;
                    break;
                }
            }
        }

        if include_boundary_gaps {
            let len_i32 = len as i32;
            let mut is = idx_start;
            let mut ie = idx_end;
            while is > 0 {
                let na = nas.get_item((is - 1) as usize)?
                    .downcast_into::<NAPosition>()?;
                if !na.borrow().is_gap {
                    break;
                }
                is -= 1;
            }
            while ie < len_i32 {
                let na = nas.get_item(ie as usize)?
                    .downcast_into::<NAPosition>()?;
                if !na.borrow().is_gap {
                    break;
                }
                ie += 1;
            }
            Ok((is as usize, ie as usize))
        } else {
            Ok((idx_start as usize, idx_end as usize))
        }
    }
}

// Private helpers (not exposed to Python).
impl NAPosition {
    fn _min_pos_list(nas: &Bound<'_, PyList>) -> PyResult<i32> {
        for i in 0..nas.len() {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            let r = na.borrow();
            if r.pos > 0 {
                return Ok(r.pos);
            }
        }
        Ok(-1)
    }

    fn _max_pos_list(nas: &Bound<'_, PyList>) -> PyResult<i32> {
        let len = nas.len();
        for i in (0..len).rev() {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            let r = na.borrow();
            if r.pos > 0 {
                return Ok(r.pos);
            }
        }
        Ok(-1)
    }

    /// Find index of element with given pos.  `first=true` → first match,
    /// `first=false` → last match.
    fn _pos2index_list(nas: &Bound<'_, PyList>, pos: i32, first: bool) -> PyResult<i32> {
        let len = nas.len();
        if first {
            for i in 0..len {
                let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
                if na.borrow().pos == pos {
                    return Ok(i as i32);
                }
            }
        } else {
            for i in (0..len).rev() {
                let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
                if na.borrow().pos == pos {
                    return Ok(i as i32);
                }
            }
        }
        Ok(-1)
    }
}

// ---------------------------------------------------------------------------
// NAPositionList — Struct-of-Arrays container
// ---------------------------------------------------------------------------

/// Iterator over `NAPositionList`, yielding `NAPosition` objects.
#[pyclass(module = "postalign_rs")]
pub struct NAPositionListIter {
    notations: Vec<u8>,
    positions: Vec<i32>,
    flags: Vec<i32>,
    index: usize,
}

#[pymethods]
impl NAPositionListIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<NAPosition> {
        if self.index >= self.notations.len() {
            return None;
        }
        let i = self.index;
        self.index += 1;
        Some(NAPosition {
            notation: self.notations[i] as i32,
            pos: self.positions[i],
            flag: self.flags[i],
            is_gap: is_gap_byte(self.notations[i]),
            payload: None,
        })
    }
}

/// Struct-of-Arrays container for nucleotide positions.
///
/// Stores `notations`, `positions`, and `flags` as contiguous Rust
/// `Vec`s, enabling vectorized bulk operations without per-element
/// Python object overhead.
#[pyclass(module = "postalign_rs")]
#[derive(Clone, Debug)]
pub struct NAPositionList {
    pub notations: Vec<u8>,
    pub positions: Vec<i32>,
    pub flags: Vec<i32>,
}

#[pymethods]
impl NAPositionList {
    // -- Core dunder methods --------------------------------------------------

    #[new]
    fn new() -> Self {
        NAPositionList {
            notations: Vec::new(),
            positions: Vec::new(),
            flags: Vec::new(),
        }
    }

    fn __len__(&self) -> usize {
        self.notations.len()
    }

    fn __bool__(&self) -> bool {
        !self.notations.is_empty()
    }

    fn __repr__(&self) -> String {
        let len = self.notations.len();
        if len <= 10 {
            let s: String = self.notations.iter().map(|&b| b as char).collect();
            format!("NAPositionList({len}, '{s}')")
        } else {
            let prefix: String = self.notations[..5].iter().map(|&b| b as char).collect();
            let suffix: String =
                self.notations[len - 5..].iter().map(|&b| b as char).collect();
            format!("NAPositionList({len}, '{prefix}...{suffix}')")
        }
    }

    fn __getitem__(&self, py: Python<'_>, index: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if let Ok(mut idx) = index.extract::<isize>() {
            let len = self.notations.len() as isize;
            if idx < 0 {
                idx += len;
            }
            if idx < 0 || idx >= len {
                return Err(pyo3::exceptions::PyIndexError::new_err(
                    "index out of range",
                ));
            }
            let i = idx as usize;
            let na = NAPosition {
                notation: self.notations[i] as i32,
                pos: self.positions[i],
                flag: self.flags[i],
                is_gap: is_gap_byte(self.notations[i]),
                payload: None,
            };
            Ok(na.into_pyobject(py)?.into_any().unbind())
        } else if let Ok(slice) = index.downcast::<PySlice>() {
            let len = self.notations.len();
            let indices = slice.indices(len as isize)?;
            let mut notations = Vec::new();
            let mut positions = Vec::new();
            let mut flags = Vec::new();
            let mut i = indices.start;
            while (indices.step > 0 && i < indices.stop)
                || (indices.step < 0 && i > indices.stop)
            {
                let idx = i as usize;
                notations.push(self.notations[idx]);
                positions.push(self.positions[idx]);
                flags.push(self.flags[idx]);
                i += indices.step;
            }
            let result = NAPositionList {
                notations,
                positions,
                flags,
            };
            Ok(result.into_pyobject(py)?.into_any().unbind())
        } else {
            Err(pyo3::exceptions::PyTypeError::new_err(
                "indices must be integers or slices",
            ))
        }
    }

    fn __iter__(&self) -> NAPositionListIter {
        NAPositionListIter {
            notations: self.notations.clone(),
            positions: self.positions.clone(),
            flags: self.flags.clone(),
            index: 0,
        }
    }

    fn __add__(&self, other: &NAPositionList) -> NAPositionList {
        let mut notations =
            Vec::with_capacity(self.notations.len() + other.notations.len());
        let mut positions =
            Vec::with_capacity(self.positions.len() + other.positions.len());
        let mut flags = Vec::with_capacity(self.flags.len() + other.flags.len());
        notations.extend_from_slice(&self.notations);
        notations.extend_from_slice(&other.notations);
        positions.extend_from_slice(&self.positions);
        positions.extend_from_slice(&other.positions);
        flags.extend_from_slice(&self.flags);
        flags.extend_from_slice(&other.flags);
        NAPositionList {
            notations,
            positions,
            flags,
        }
    }

    // -- Constructors --------------------------------------------------------

    /// Create from a byte sequence (replaces init_from_bytes).
    #[classmethod]
    fn from_bytes(
        _cls: &Bound<'_, pyo3::types::PyType>,
        seq_text: &[u8],
    ) -> NAPositionList {
        let len = seq_text.len();
        let mut notations = Vec::with_capacity(len);
        let mut positions = Vec::with_capacity(len);
        let flags = vec![0i32; len];
        let mut offset: i32 = 1;
        for &b in seq_text {
            let upper = b.to_ascii_uppercase();
            notations.push(upper);
            if is_gap_byte(upper) {
                positions.push(-1);
            } else {
                positions.push(offset);
                offset += 1;
            }
        }
        NAPositionList {
            notations,
            positions,
            flags,
        }
    }

    /// Create from flat arrays.
    #[classmethod]
    fn from_arrays(
        _cls: &Bound<'_, pyo3::types::PyType>,
        notations: Vec<u8>,
        positions: Vec<i32>,
        flags: Vec<i32>,
    ) -> PyResult<NAPositionList> {
        let len = notations.len();
        if positions.len() != len || flags.len() != len {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "notations, positions, and flags must have equal length",
            ));
        }
        Ok(NAPositionList {
            notations,
            positions,
            flags,
        })
    }

    /// Create a list of gap positions.
    #[classmethod]
    #[pyo3(signature = (gaplen,))]
    fn init_gaps(
        _cls: &Bound<'_, pyo3::types::PyType>,
        gaplen: usize,
    ) -> NAPositionList {
        NAPositionList {
            notations: vec![GAP_DASH; gaplen],
            positions: vec![-1; gaplen],
            flags: vec![0; gaplen],
        }
    }

    /// Convert from list[NAPosition] to NAPositionList.
    #[classmethod]
    fn from_list(
        _cls: &Bound<'_, pyo3::types::PyType>,
        nas: &Bound<'_, PyList>,
    ) -> PyResult<NAPositionList> {
        let len = nas.len();
        let mut notations = Vec::with_capacity(len);
        let mut positions = Vec::with_capacity(len);
        let mut flags = Vec::with_capacity(len);
        for i in 0..len {
            let na = nas.get_item(i)?.downcast_into::<NAPosition>()?;
            let r = na.borrow();
            notations.push(r.notation as u8);
            positions.push(r.pos);
            flags.push(r.flag);
        }
        Ok(NAPositionList {
            notations,
            positions,
            flags,
        })
    }

    /// Convert to list[NAPosition].
    fn to_list(&self) -> Vec<NAPosition> {
        (0..self.notations.len())
            .map(|i| NAPosition {
                notation: self.notations[i] as i32,
                pos: self.positions[i],
                flag: self.flags[i],
                is_gap: is_gap_byte(self.notations[i]),
                payload: None,
            })
            .collect()
    }

    // -- Bulk read operations ------------------------------------------------

    fn count_gaps(&self) -> usize {
        self.notations.iter().filter(|&&b| is_gap_byte(b)).count()
    }

    fn count_nongaps(&self) -> usize {
        self.notations
            .iter()
            .filter(|&&b| !is_gap_byte(b))
            .count()
    }

    fn any_has_gap(&self) -> bool {
        self.notations.iter().any(|&b| is_gap_byte(b))
    }

    fn all_have_gap(&self) -> bool {
        self.notations.iter().all(|&b| is_gap_byte(b))
    }

    fn any_has_flag(&self, flag: i32) -> bool {
        self.flags.iter().any(|&f| (f & flag) != 0)
    }

    fn all_have_flag(&self, flag: i32) -> bool {
        self.flags.iter().all(|&f| (f & flag) != 0)
    }

    fn as_bytes<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.notations)
    }

    fn as_str(&self) -> String {
        self.notations.iter().map(|&b| b as char).collect()
    }

    // -- Bulk write ----------------------------------------------------------

    /// Set flag (bitwise OR) on every element, in-place.
    fn set_flag(&mut self, flag: i32) {
        for f in &mut self.flags {
            *f |= flag;
        }
    }

    // -- Index operations ----------------------------------------------------

    fn min_pos(&self) -> i32 {
        for &p in &self.positions {
            if p > 0 {
                return p;
            }
        }
        -1
    }

    fn max_pos(&self) -> i32 {
        for &p in self.positions.iter().rev() {
            if p > 0 {
                return p;
            }
        }
        -1
    }

    #[pyo3(signature = (start=-1, stop=-1))]
    fn min_nongap_index(&self, start: i32, stop: i32) -> i32 {
        let len = self.positions.len() as i32;
        for idx in 0..len {
            if start > -1 && idx < start {
                continue;
            }
            if stop > -1 && idx >= stop {
                break;
            }
            if self.positions[idx as usize] > 0 {
                return idx;
            }
        }
        -1
    }

    #[pyo3(signature = (start=-1, stop=-1))]
    fn max_nongap_index(&self, start: i32, stop: i32) -> i32 {
        let len = self.positions.len() as i32;
        for revidx in 0..len {
            let idx = len - 1 - revidx;
            if start > -1 && idx < start {
                break;
            }
            if stop > -1 && idx >= stop {
                continue;
            }
            if self.positions[idx as usize] > 0 {
                return idx;
            }
        }
        -1
    }

    fn remove_gaps(&self) -> NAPositionList {
        let len = self.notations.len();
        let mut notations = Vec::with_capacity(len);
        let mut positions = Vec::with_capacity(len);
        let mut flags = Vec::with_capacity(len);
        for i in 0..len {
            if !is_gap_byte(self.notations[i]) {
                notations.push(self.notations[i]);
                positions.push(self.positions[i]);
                flags.push(self.flags[i]);
            }
        }
        NAPositionList {
            notations,
            positions,
            flags,
        }
    }

    #[pyo3(signature = (pos_start, pos_end, include_boundary_gaps=false))]
    fn posrange2indexrange(
        &self,
        pos_start: i32,
        pos_end: i32,
        include_boundary_gaps: bool,
    ) -> (usize, usize) {
        let max_p = self._max_pos_inner();
        let min_p = self._min_pos_inner();
        if max_p < 0 || min_p < 0 {
            return (0, 0);
        }

        let idx_start: i64;
        let idx_end: i64;

        if pos_start > max_p {
            let v = self._pos2index_last(max_p);
            idx_start = v;
            idx_end = v;
        } else if pos_end < min_p {
            let v = self._pos2index_first(min_p);
            idx_start = v;
            idx_end = v;
        } else {
            let ps = if min_p > pos_start { min_p } else { pos_start };
            let pe = if max_p < pos_end { max_p } else { pos_end };

            let mut is = -1i64;
            for pos in ps..=pe {
                is = self._pos2index_first(pos);
                if is > -1 {
                    break;
                }
            }
            idx_start = is;

            let mut ie = -1i64;
            for pos in (ps..=pe).rev() {
                let v = self._pos2index_last(pos);
                if v > -1 {
                    ie = v + 1;
                    break;
                }
            }
            idx_end = ie;
        }

        if include_boundary_gaps {
            let n = self.notations.len() as i64;
            let mut is = idx_start;
            let mut ie = idx_end;
            while is > 0 && is_gap_byte(self.notations[(is - 1) as usize]) {
                is -= 1;
            }
            while ie < n && is_gap_byte(self.notations[ie as usize]) {
                ie += 1;
            }
            (is as usize, ie as usize)
        } else {
            (idx_start as usize, idx_end as usize)
        }
    }
}

// Private helpers (not exposed to Python).
impl NAPositionList {
    fn _min_pos_inner(&self) -> i32 {
        for &p in &self.positions {
            if p > 0 {
                return p;
            }
        }
        -1
    }

    fn _max_pos_inner(&self) -> i32 {
        for &p in self.positions.iter().rev() {
            if p > 0 {
                return p;
            }
        }
        -1
    }

    fn _pos2index_first(&self, pos: i32) -> i64 {
        for (idx, &p) in self.positions.iter().enumerate() {
            if p == pos {
                return idx as i64;
            }
        }
        -1
    }

    fn _pos2index_last(&self, pos: i32) -> i64 {
        for idx in (0..self.positions.len()).rev() {
            if self.positions[idx] == pos {
                return idx as i64;
            }
        }
        -1
    }
}

// ---------------------------------------------------------------------------
// AAPosition (stub — all methods raise NotImplementedError)
// ---------------------------------------------------------------------------

/// Amino-acid position stub.
#[pyclass(module = "postalign_rs")]
#[derive(Debug)]
pub struct AAPosition {
    #[pyo3(get, set)]
    pub notation: i32,
    #[pyo3(get, set)]
    pub pos: i32,
    #[pyo3(get, set)]
    pub flag: i32,
    #[pyo3(get)]
    pub is_gap: bool,
    #[pyo3(get, set)]
    pub payload: Option<PyObject>,
}

impl Clone for AAPosition {
    fn clone(&self) -> Self {
        AAPosition {
            notation: self.notation,
            pos: self.pos,
            flag: self.flag,
            is_gap: self.is_gap,
            payload: clone_opt_pyobj(&self.payload),
        }
    }
}

#[pymethods]
impl AAPosition {
    #[new]
    #[pyo3(signature = (notation, pos, flag, payload=None))]
    fn new(notation: i32, pos: i32, flag: i32, payload: Option<PyObject>) -> Self {
        AAPosition {
            notation,
            pos,
            flag,
            is_gap: is_gap_byte(notation as u8),
            payload,
        }
    }

    fn __copy__(&self) -> PyResult<Self> {
        Err(pyo3::exceptions::PyNotImplementedError::new_err(
            "Amino acid sequence is not yet supported",
        ))
    }

    #[classmethod]
    #[pyo3(signature = (gaplen,))]
    fn init_gaps(
        _cls: &Bound<'_, pyo3::types::PyType>,
        gaplen: usize,
    ) -> PyResult<Vec<AAPosition>> {
        let _ = gaplen;
        Err(pyo3::exceptions::PyNotImplementedError::new_err(
            "Amino acid sequence is not yet supported",
        ))
    }

    #[classmethod]
    #[pyo3(signature = (seq_text, seq_payload=None))]
    fn init_from_bytes(
        _cls: &Bound<'_, pyo3::types::PyType>,
        seq_text: &[u8],
        seq_payload: Option<Vec<PyObject>>,
    ) -> PyResult<Vec<AAPosition>> {
        let _ = (seq_text, seq_payload);
        Err(pyo3::exceptions::PyNotImplementedError::new_err(
            "Amino acid sequence is not yet supported",
        ))
    }
}

// ---------------------------------------------------------------------------
// Standalone function: enumerate_seq_pos
// ---------------------------------------------------------------------------

/// Assign sequential 1-based positions to non-gap bytes, -1 for gaps.
#[pyfunction]
pub fn enumerate_seq_pos(seq_text: &[u8]) -> Vec<i32> {
    let mut offset: i32 = 1;
    seq_text
        .iter()
        .map(|&b| {
            if is_gap_byte(b) {
                -1
            } else {
                let p = offset;
                offset += 1;
                p
            }
        })
        .collect()
}
