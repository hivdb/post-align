# post-align

A post alignment toolkit for refining pairwise/multiple alignment sequences.

## Features

### Codon alignment

Realigns gaps within coding regions to respect codon boundaries, improving
downstream translation and mutation calling accuracy.

### Adjust gap placement

Repositions gaps to maximize BLOSUM62 amino-acid similarity scores, using
configurable gap placement scoring.

### Generate paired codon / mutation report

## Performance

The codon alignment engine has been progressively optimized across five tiers:

| Tier | Implementation | Throughput | Speedup vs T1 |
|------|----------------|------------|----------------|
| T1 | Python (original) | ~126 reads/sec | 1.0x |
| T2 | Python (optimized) | ~162 reads/sec | 1.3x |
| T3 | Cython (byte-array) | ~198 reads/sec | 1.6x |
| T4 | Rust/PyO3 | ~728 reads/sec | 5.8x |
| T5 | Rust optimized (single) | ~944 reads/sec | 7.5x |
| T5 | Rust batch + rayon (8 cores) | ~5,081 reads/sec | ~5x vs T1×8¹ |

¹ Fair comparison: T5 batch uses 8 cores via rayon. Comparing against
T1 × 8 cores (theoretical ~1,008 reads/sec), the per-core speedup is ~5x.

T5 optimizations include:
- Precomputed 128×128 IUPAC nucleotide scoring table
- Inline codon translation (zero heap allocation)
- Incremental per-codon score caching
- Pre-allocated SoA work buffers
- Full boundary detection in Rust (Phase C)
- Batch API with rayon parallelism (Phase D)

## Development

### Prerequisites

- Python >= 3.13
- Rust toolchain (for T4/T5)
- pipenv

### Setup

```bash
pipenv install --dev
```

### Build extensions

```bash
# Cython extensions (T2, T3)
make build-ext

# Rust extension (T4, T5)
pipenv run maturin develop --release --manifest-path rust_ext/Cargo.toml
```

### Run tests

```bash
make test
```

### Run benchmarks

```bash
make benchmark
```
