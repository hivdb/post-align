# Codon Alignment — Overview

## Purpose

After an initial pairwise alignment (e.g. NW / SW), gap positions are
often arbitrary with respect to the reading frame.  The **codon
alignment** post-processor slides every gap block to the position that
maximises a combined IUPAC nucleotide + BLOSUM62 amino-acid score,
ensuring that insertions and deletions land on codon boundaries whenever
possible.

This is critical for HIV drug-resistance genotyping pipelines where a
single misplaced gap can shift the reading frame and produce incorrect
amino-acid calls downstream.

## High-Level Pipeline

```
input pair (ref, seq)
        │
        ▼
  ┌──────────────┐
  │  Boundary    │  Clip to ref_start..ref_end; align to codon
  │  Detection   │  boundary (pos % 3 == 0 relative to ref_start).
  └──────┬───────┘
         ▼
  ┌──────────────┐
  │  Gather Gaps │  Merge gap windows within min_gap_distance of each
  │              │  other.  Remove redundant paired gaps.  Centre
  │              │  remaining gaps within each window.
  └──────┬───────┘
         ▼
  ┌──────────────┐
  │  Group by    │  Split the paired sequences into codon-aligned
  │  Codons      │  groups (new group every 3 non-gap ref bases).
  └──────┬───────┘
         ▼
  ┌──────────────┐
  │  Adjust Gap  │  For each contiguous run of gap-containing codons:
  │  Placement   │    1. Extend a flanking window (±window_size codons).
  │              │    2. Try every codon-aligned position (centre-expand
  │              │       search from the original gap location).
  │              │    3. Score each candidate with IUPAC + BLOSUM62.
  │              │    4. Apply gap_placement_score bonuses.
  │              │    5. Break ties: codon-boundary > original > leftmost.
  │              │    6. Keep the best.
  └──────┬───────┘
         ▼
  ┌──────────────┐
  │  Finalise    │  Move gaps to codon ends (cosmetic) and return the
  │              │  updated (ref, seq) Sequence objects.
  └──────────────┘
```

## Backends

| Backend | Module | When to use |
|---------|--------|-------------|
| **Python/Cython** | `postalign.processors.codon_alignment` | Default for development, debugging, and when the Rust extension is unavailable. Compiled to C via Cython for ~2× speedup over pure Python. |
| **Rust** | `postalign.processors.codon_alignment_rust` | Production workloads (NGS scale).  Same algorithm implemented entirely in Rust via PyO3 with additional optimisations (compile-time IUPAC table, stack-allocated amino-acid translation, incremental per-codon score cache). |

Select the backend at the CLI level:

```bash
postalign codon-alignment --backend rust   # default
postalign codon-alignment --backend python
```

Or programmatically:

```python
from postalign.processors.codon_alignment import codon_align       # Python
from postalign.processors.codon_alignment_rust import codon_align   # Rust
```

Both expose the identical `codon_align(refseq, seq, ...)` signature and
produce identical outputs for the same inputs.

## Scoring

Each candidate gap position is scored as:

```
score = base_gap_penalty
      + Σ  iupac_score(my[i], other[i])      (per nucleotide)
      + Σ  blosum62_score(my_aa[c], other_aa[c])  (per codon)
      + gap_placement_bonus(pos, gap_len)     (optional)
```

- **Base gap penalty**: `−gap_length`, zeroed when the gap is at the
  sequence start or end.
- **IUPAC score**: `1.0` for identity, fractional for partial IUPAC
  overlap, negative for mismatch.
- **BLOSUM62 score**: standard amino-acid substitution score, averaged
  over ambiguous translations.
- **Gap placement bonus**: user-supplied `{(position, gap_len): score}`
  overrides parsed from the `--gap-placement-score` CLI option.

## Key Parameters

| Parameter | CLI flag | Default | Description |
|-----------|----------|---------|-------------|
| `min_gap_distance` | `--min-gap-distance` | 30 | Merge gap windows within this many NA positions. |
| `window_size` | `--window-size` | 10 | Flanking codons used when searching for the best gap position. |
| `gap_placement_score` | `--gap-placement-score` | `{}` | Per-position scoring overrides (see CLI `--help`). |
| `ref_start` / `ref_end` | `--ref-start` / `--ref-end` | full range | Reference position range to process. |

## Performance

The Rust backend achieves **30–80× speedup** over pure Python on
typical HIV sequences (~3 kb), and the batch API
(`postalign_rs.codon_align_batch`) adds rayon-based parallelism for
multi-read workloads.
