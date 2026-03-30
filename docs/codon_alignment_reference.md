# Codon Alignment — Technical Reference

## Module Layout

```
postalign/processors/
├── codon_alignment.py          # Python/Cython backend + CLI entry point
├── codon_alignment_rust.py     # Rust backend (thin PyO3 wrapper)
└── __init__.py                 # re-exports codon_alignment

rust_ext/src/
├── lib.rs                      # PyO3 bindings, batch API, boundary detection
├── align.rs                    # Core alignment algorithm (gap gathering,
│                               #   codon grouping, gap placement search)
└── scoring.rs                  # IUPAC table, BLOSUM62 table, codon translation,
                                #   InlineAA, scoring functions
```

---

## Python API

### `codon_alignment.codon_align`

```python
def codon_align(
    refseq: Sequence,
    seq: Sequence,
    min_gap_distance: int,
    window_size: int,
    gap_placement_score: dict[int, dict[tuple[int, int], int]],
    ref_start: int,
    ref_end: int,
) -> RefSeqPair:
```

Main entry point for the **Python/Cython** backend.

| Parameter | Type | Description |
|-----------|------|-------------|
| `refseq` | `Sequence` | Reference sequence (immutable `NAPosition` list). |
| `seq` | `Sequence` | Target sequence to realign. |
| `min_gap_distance` | `int` | Merge gap windows within this NA distance. |
| `window_size` | `int` | Flanking codons for gap search. |
| `gap_placement_score` | `dict` | `{gap_type: {(pos, size): score}}`. |
| `ref_start` | `int` | First reference position (1-based, inclusive). |
| `ref_end` | `int` | Last reference position (1-based, inclusive). |

**Returns** `(refseq, seq)` — updated `Sequence` objects with a new
`seqtext` layer tagged `codonalign(ref_start,ref_end)`.

### `codon_alignment_rust.codon_align`

Identical signature and return type.  Internally marshals `NAPosition`
lists to flat arrays, calls `postalign_rs.codon_align_full`, and
reconstructs `Sequence` objects from the result.

---

## Internal Functions (Python backend)

### Gap Gathering

| Function | Description |
|----------|-------------|
| `_find_windows_with_gap(refnas, seqnas, min_gap_distance)` | Scan paired NAs; merge gap indices within `min_gap_distance` into `(start, end)` windows. |
| `_remove_redundant_gaps(refnas, seqnas)` | Strip `min(ref_gaps, seq_gaps)` paired gaps. |
| `_move_gaps_to_center(nas)` | Separate gaps from non-gaps, re-insert gaps at centre. |
| `_gather_gaps(refnas, seqnas, min_gap_distance)` | Orchestrates the above per window, right-to-left. |

### Codon Grouping

| Function | Description |
|----------|-------------|
| `group_by_codons(refnas, seqnas)` | Split into codon-aligned `(refcodons, seqcodons)` lists. New group every 3 non-gap ref bases. |
| `_find_codon_trim_slice(codons)` | Trim leading/trailing all-gap codons (boundary gaps). |

### Gap Placement Search

| Function | Description |
|----------|-------------|
| `_find_best_matches(mynas, othernas, bp1_indices, gap_type, gps, is_start, is_end)` | Try every codon-aligned position for the gap block using centre-expand order.  Compute IUPAC + BLOSUM62 score at each, apply GPS bonuses, select the best with tie-breaking. |
| `_paired_find_best_matches(refnas, seqnas, gap_type, gps, is_start, is_end)` | Compute `bp1_indices` then dispatch to `_find_best_matches` for the gapped side. |
| `_adjust_gap_placement(refcodons, seqcodons, window_size, gps, is_start, is_end)` | Iterate gap groups, extend windows, call `_paired_find_best_matches`, splice results back. |

### Top-Level

| Function | Description |
|----------|-------------|
| `realign_gaps(refnas, seqnas, min_gap_distance, window_size, gps, is_start, is_end)` | `gather_gaps` → `group_by_codons` → `adjust_gap_placement` → `move_gap_to_codon_end` → flatten. |
| `codon_align(...)` | Boundary detection → `realign_gaps` → `push_seqtext`. |

---

## Rust API (`postalign_rs`)

All functions are exposed as a flat PyO3 module.

### `postalign_rs.realign_gaps`

Baseline Rust port of `realign_gaps`.  Same flat-array interface.

### `postalign_rs.realign_gaps_optimized`

Optimised Rust `realign_gaps` with:

- **Compile-time IUPAC table** (`scoring::IUPAC_TABLE`, 128×128 `i32`
  array built at `const` time).
- **`InlineAA`** — stack-allocated `[u8; 4]` codon translation (no heap
  per codon).
- **Incremental per-codon score cache** — only recompute codons in the
  changed range when the gap slides by one step.
- **Centre-expand search** — start from the original gap position and
  spiral outward, enabling early-exit heuristics in future.

### `postalign_rs.codon_align_full`

Full pipeline in Rust: boundary detection + `realign_gaps_optimized`.
Returns `None` or `(idx_start, idx_end, ref_n, ref_p, ref_f, seq_n,
seq_p, seq_f)`.

### `postalign_rs.codon_align_batch`

Batch API: accepts a Python list of 8-tuples, releases the GIL, and
processes all pairs in parallel via **rayon**.  Returns a list of
results matching `codon_align_full` output.

---

## Rust Internals

### `align.rs`

| Item | Description |
|------|-------------|
| `NaPos` | `{notation: u8, pos: i32, flag: i32}` — mirrors `NAPosition`. |
| `GapPlacementScore` | `{refgap, seqgap}` hash maps keyed by `(pos, gap_len)`. |
| `gather_gaps` | Merge windows, remove redundant gaps, centre remaining gaps. |
| `group_by_codons` | Split into codon groups by reference reading frame. |
| `adjust_gap_placement` / `adjust_gap_placement_optimized` | Gap search using baseline / optimised scoring. |
| `realign_gaps` / `realign_gaps_optimized` | Top-level orchestrators. |
| `find_best_matches` / `find_best_matches_optimized` | Per-gap-group scoring loop (baseline / optimised). |
| `center_expand_positions` | Generate search order spiralling outward from the original gap index. |

### `scoring.rs`

| Item | Description |
|------|-------------|
| `IUPAC_TABLE` | `[i32; 16384]` — compile-time IUPAC score × 1000.  Lookup: `iupac_score(a, b)`. |
| `BLOSUM62_TABLE` | `[i8; 16384]` — compile-time BLOSUM62 lookup by ASCII byte. |
| `InlineAA` | `{data: [u8; 4], len: u8}` — stack amino-acid set. |
| `translate_codon_inline` | IUPAC-aware codon → `InlineAA` (handles gaps, frameshifts, ambiguity). |
| `compute_full_score` | Full-window IUPAC + BLOSUM62 with per-codon vectors for caching. |
| `recompute_codon_score` | Single-codon delta recomputation for incremental updates. |
| `blosum62_score` / `blosum62_score_inline` | AA-pair scoring (Vec / InlineAA variants). |

### `lib.rs`

| Item | Description |
|------|-------------|
| `parse_gap_placement_score` | Python dict → `GapPlacementScore`. |
| `posrange2indexrange` | Boundary detection with gap inclusion. |
| `codon_align_core` | Pure-Rust pipeline (no PyO3 types) used by both single and batch APIs. |
| `BatchInput` | Owned struct for rayon `Send + Sync` safety. |

---

## CLI

The `codon-alignment` command is registered via Click in
`codon_alignment.py`:

```
postalign codon-alignment [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--min-gap-distance` | `INT` | 30 | Merge gap windows within this distance. |
| `--window-size` | `INT` | 10 | Flanking codons for search. |
| `--gap-placement-score` | `TEXT` | `""` | Comma-separated `pos/sizeType:score` tokens. |
| `--backend` | `python\|rust` | `rust` | Execution backend. |

The command returns a `Processor` callable that is applied to each
`(refseq, seq)` pair in the post-alignment pipeline.

---

## Constants

| Name | Value | Meaning |
|------|-------|---------|
| `REFGAP` | `1` | Gap type: gap in reference. |
| `SEQGAP` | `2` | Gap type: gap in sequence. |
| `NOGAP` | `0` | No gap in either side. |

---

## Testing

```bash
# Unit + property-based tests (both backends)
pipenv run pytest tests/test_codon_alignment.py tests/test_hypothesis.py -v

# With coverage
pipenv run pytest tests/ --cov=postalign.processors.codon_alignment \
                         --cov=postalign.processors.codon_alignment_rust \
                         --cov-report=term-missing
```

Tests are parametrised over both `python` and `rust` backends via an
autouse fixture.  Property-based tests use Hypothesis to verify
structural invariants (length preservation, non-gap content
preservation) on random sequences.
