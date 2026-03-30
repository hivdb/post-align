# NAPositionList Refactor Plan

**Goal:** Replace `list[NAPosition]` with a Struct-of-Arrays `NAPositionList` for better performance.

**Expected Gains:**
- Eliminate per-element Python object overhead (10K objects → 1 object)
- Zero-copy data passing to Rust (no marshaling)
- Vectorized bulk operations (~50-2000× faster for set_flag, count_gaps, etc.)

---

## Phase 1: Create NAPositionList in Rust

**Files:** `rust_ext/src/position.rs`, `rust_ext/src/lib.rs`

- [ ] **1.0** Add `NAPositionList` struct with SoA storage
- [ ] **1a** Implement core methods: `__len__`, `__getitem__`, `__iter__`, `from_bytes`, `from_arrays`
- [ ] **1b** Implement bulk ops: `count_gaps`, `count_nongaps`, `any_has_gap`, `all_have_gap`, `set_flag`, `as_bytes`, `as_str`
- [ ] **1c** Implement index ops: `min_pos`, `max_pos`, `min_nongap_index`, `max_nongap_index`, `posrange2indexrange`, `remove_gaps`
- [ ] **1d** Add `init_gaps` classmethod and slicing support

---

## Phase 2: Integrate with Rust Codon Alignment

**Files:** `rust_ext/src/lib.rs`, `postalign/processors/codon_alignment_rust.py`

- [ ] **2.0** Update `codon_alignment_rust.py` to use `NAPositionList`
- [ ] **2a** Modify `codon_align_full` to accept `NAPositionList` directly (zero-copy)
- [ ] **2b** Return `NAPositionList` from Rust (eliminate reconstruction loop)

---

## Phase 3: Update Parsers

**Files:** `postalign/parsers/fasta.py`, `postalign/parsers/paf.py`

- [ ] **3a** Update `fasta.py`: `init_from_bytes` → `NAPositionList.from_bytes`
- [ ] **3b** Update `paf.py`: use `NAPositionList` for reftext/seqtext

---

## Phase 4: Update Sequence Model & Consumers

**Files:** `postalign/models/sequence.py`, `postalign/models/_sequence.py`, `postalign/processors/*.py`, `postalign/utils/*.py`

- [ ] **4a** Change `Sequence.seqtext` type from `list[Position]` to `NAPositionList`
- [ ] **4b** Update `codon_alignment.py` (pure Python path)
- [ ] **4c** Update `save_json.py`, `save_fasta.py`, `trim_by_ref.py`
- [ ] **4d** Update `group_by_codons.py` and other utils

---

## Phase 5: Testing & Validation

- [ ] **5.0** Run full test suite (385 tests), fix regressions
- [ ] **5.1** Add NAPositionList-specific tests

---

## Phase 6: Benchmarking

- [ ] **6.0** Run `benchmarks/compare_cython_rust.py` before/after
- [ ] **6.1** Document performance improvements

---

## Phase 7: Cleanup

- [ ] **7.0** Remove unused `payload` field from NAPosition
- [ ] **7.1** Deprecate old `list[NAPosition]` patterns
- [ ] **7.2** Update type hints and documentation

---

## Notes

- `NAPositionList` can coexist with `list[NAPosition]` during migration via `__iter__`
- **Phase 2 provides the biggest performance win** (zero-copy Rust calls)
- Phases 3-4 can be done incrementally

---

## Current Status

**Started:** Not yet  
**Last Updated:** $(date +%Y-%m-%d)
