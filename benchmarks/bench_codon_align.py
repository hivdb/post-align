"""Benchmark codon_align() across all five tiers (T1–T5).

Run with:  pipenv run python -m pytest benchmarks/ -v --benchmark-only
"""
import pytest
import postalign_rs

from postalign.processors.codon_alignment import (
    codon_align as t1_codon_align,
)
from postalign.processors.codon_alignment_optimized import (
    codon_align as t2_codon_align,
)
from postalign.processors.codon_alignment_cython import (
    codon_align as t3_codon_align,
)
from postalign.processors.codon_alignment_rust import (
    codon_align as t4_codon_align,
)
from postalign.processors.codon_alignment_rust_t5 import (
    codon_align as t5_codon_align,
)

from benchmarks.conftest import make_pair, DEFAULT_GPS, SCENARIOS
from postalign.models.na_position import NAPosition


_IMPLS = {
    'T1_original': t1_codon_align,
    'T2_optim':    t2_codon_align,
    'T3_cython':   t3_codon_align,
    'T4_rust':     t4_codon_align,
    'T5_rust_opt': t5_codon_align,
}


def _make_args(length, n_ref_gaps, n_seq_gaps, gap_size):
    ref, seq, gps = make_pair(
        length, n_ref_gaps, n_seq_gaps, gap_size, seed=42)
    ref_start = 1
    ref_end = NAPosition.max_pos(ref.seqtext)
    return ref, seq, gps, ref_start, ref_end


# ---------------------------------------------------------------------------
# Single-read latency benchmarks — one per (tier × scenario)
# ---------------------------------------------------------------------------

def _bench_ids():
    ids = []
    for sc in SCENARIOS:
        for impl in _IMPLS:
            ids.append(f'{impl}-{sc}')
    return ids


def _bench_params():
    params = []
    for sc_name, sc_kw in SCENARIOS.items():
        for impl_name, impl_fn in _IMPLS.items():
            params.append((impl_name, impl_fn, sc_name, sc_kw))
    return params


@pytest.mark.parametrize(
    'impl_name,impl_fn,sc_name,sc_kw',
    _bench_params(),
    ids=_bench_ids(),
)
def test_single_read_latency(
    benchmark, impl_name, impl_fn, sc_name, sc_kw
):
    ref, seq, gps, ref_start, ref_end = _make_args(**sc_kw)
    benchmark.group = sc_name
    benchmark.extra_info['tier'] = impl_name
    benchmark.extra_info['scenario'] = sc_name

    benchmark(
        impl_fn,
        ref, seq, 30, 10, gps, ref_start, ref_end,
    )


# ---------------------------------------------------------------------------
# Throughput benchmarks — N reads through each tier
# ---------------------------------------------------------------------------

THROUGHPUT_N = 500


def _throughput_pairs():
    """Pre-generate N pairs for throughput measurement."""
    pairs = []
    kw = SCENARIOS['medium_multi']
    for i in range(THROUGHPUT_N):
        ref, seq, gps = make_pair(**kw, seed=i)
        ref_start = 1
        ref_end = NAPosition.max_pos(ref.seqtext)
        pairs.append((ref, seq, gps, ref_start, ref_end))
    return pairs


_PAIRS = _throughput_pairs()


def _run_batch(impl_fn):
    for ref, seq, gps, rs, re in _PAIRS:
        impl_fn(ref, seq, 30, 10, gps, rs, re)


@pytest.mark.parametrize(
    'impl_name,impl_fn',
    list(_IMPLS.items()),
    ids=list(_IMPLS.keys()),
)
def test_throughput_500(benchmark, impl_name, impl_fn):
    benchmark.group = 'throughput_500'
    benchmark.extra_info['tier'] = impl_name
    benchmark.extra_info['n_reads'] = THROUGHPUT_N
    benchmark(lambda: _run_batch(impl_fn))


# ---------------------------------------------------------------------------
# Batch throughput benchmark — rayon parallel via codon_align_batch
# ---------------------------------------------------------------------------


def _prepare_batch_items():
    """Pre-extract flat arrays for batch Rust API."""
    items = []
    for ref, seq, gps, rs, re in _PAIRS:
        rn = [na.notation for na in ref.seqtext]
        rp = [na.pos for na in ref.seqtext]
        rf = [na.flag for na in ref.seqtext]
        sn = [na.notation for na in seq.seqtext]
        sp = [na.pos for na in seq.seqtext]
        sf = [na.flag for na in seq.seqtext]
        items.append((rn, rp, rf, sn, sp, sf, rs, re))
    return items, gps


_BATCH_ITEMS, _BATCH_GPS = _prepare_batch_items()


def test_throughput_500_batch_rayon(benchmark):
    benchmark.group = 'throughput_500'
    benchmark.extra_info['tier'] = 'T5_batch_rayon'
    benchmark.extra_info['n_reads'] = THROUGHPUT_N
    benchmark(
        postalign_rs.codon_align_batch,
        _BATCH_ITEMS, 30, 10, _BATCH_GPS,
    )
