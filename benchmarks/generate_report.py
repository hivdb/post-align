#!/usr/bin/env python
"""Generate benchmark comparison charts from pytest-benchmark JSON output.

Usage:
    # 1) Run benchmarks and save JSON:
    pipenv run python -m pytest benchmarks/ --benchmark-only \
        --benchmark-json=benchmarks/data/results.json

    # 2) Generate graphs:
    pipenv run python benchmarks/generate_report.py
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402


DATA_DIR = Path(__file__).parent / 'data'
RESULTS_FILE = DATA_DIR / 'results.json'
OUTPUT_DIR = DATA_DIR / 'charts'

TIER_ORDER = [
    'T1_original', 'T2_optim', 'T3_cython',
    'T4_rust', 'T5_rust_opt', 'T5_batch_rayon',
]
TIER_LABELS = {
    'T1_original': 'T1 Original',
    'T2_optim': 'T2 Optimized',
    'T3_cython': 'T3 Cython',
    'T4_rust': 'T4 Rust',
    'T5_rust_opt': 'T5 Rust Opt',
    'T5_batch_rayon': 'T5 Batch+Rayon',
}
TIER_COLORS = {
    'T1_original': '#4e79a7',
    'T2_optim': '#f28e2b',
    'T3_cython': '#59a14f',
    'T4_rust': '#e15759',
    'T5_rust_opt': '#76b7b2',
    'T5_batch_rayon': '#b07aa1',
}
# Tiers that participate in latency charts (batch excluded)
LATENCY_TIERS = [
    'T1_original', 'T2_optim', 'T3_cython',
    'T4_rust', 'T5_rust_opt',
]


def load_results(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def parse_benchmarks(data: dict) -> dict:
    """Parse pytest-benchmark JSON into structured results."""
    results = defaultdict(lambda: defaultdict(dict))
    for bench in data['benchmarks']:
        name = bench['name']
        extra = bench.get('extra_info', {})
        tier = extra.get('tier', '')
        scenario = extra.get('scenario', '')
        n_reads = extra.get('n_reads', 0)
        stats = bench['stats']

        if 'throughput' in name:
            results['throughput'][tier] = {
                'median': stats['median'],
                'mean': stats['mean'],
                'stddev': stats['stddev'],
                'n_reads': n_reads,
                'rounds': stats['rounds'],
            }
        elif scenario:
            results['latency'][scenario][tier] = {
                'median': stats['median'],
                'mean': stats['mean'],
                'min': stats['min'],
                'max': stats['max'],
                'stddev': stats['stddev'],
                'rounds': stats['rounds'],
            }
    return dict(results)


def plot_latency_bars(results: dict, output_dir: Path):
    """Grouped bar chart: latency per tier at each scenario."""
    scenarios = list(results['latency'].keys())
    if not scenarios:
        return

    fig, ax = plt.subplots(figsize=(12, 6))
    n_scenarios = len(scenarios)
    n_tiers = len(LATENCY_TIERS)
    bar_width = 0.15
    x_base = range(n_scenarios)

    for i, tier in enumerate(LATENCY_TIERS):
        medians = []
        for sc in scenarios:
            val = results['latency'].get(sc, {}).get(tier, {})
            medians.append(val.get('median', 0) * 1000)  # ms
        offset = (i - n_tiers / 2 + 0.5) * bar_width
        bars = ax.bar(
            [x + offset for x in x_base], medians,
            bar_width, label=TIER_LABELS.get(tier, tier),
            color=TIER_COLORS.get(tier, '#999'))
        for bar, val in zip(bars, medians):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                        f'{val:.2f}', ha='center', va='bottom', fontsize=7)

    ax.set_xlabel('Scenario')
    ax.set_ylabel('Median latency (ms)')
    ax.set_title('Single-Read Latency by Scenario and Tier')
    ax.set_xticks(list(x_base))
    ax.set_xticklabels(scenarios, rotation=15)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_dir / 'latency_bars.png', dpi=150)
    plt.close(fig)
    print(f'  -> {output_dir / "latency_bars.png"}')


def plot_speedup(results: dict, output_dir: Path):
    """Speedup chart: T2/T3/T4 relative to T1 baseline."""
    scenarios = list(results['latency'].keys())
    if not scenarios:
        return

    fig, ax = plt.subplots(figsize=(12, 6))
    n_scenarios = len(scenarios)
    bar_width = 0.22
    x_base = range(n_scenarios)

    for i, tier in enumerate(LATENCY_TIERS[1:], 0):  # skip T1
        speedups = []
        for sc in scenarios:
            t1_val = results['latency'].get(sc, {}).get(
                'T1_original', {}).get('median', 1)
            tier_val = results['latency'].get(sc, {}).get(
                tier, {}).get('median', 1)
            speedups.append(t1_val / tier_val if tier_val > 0 else 0)
        offset = (i - 1) * bar_width
        bars = ax.bar(
            [x + offset for x in x_base], speedups,
            bar_width, label=TIER_LABELS.get(tier, tier),
            color=TIER_COLORS.get(tier, '#999'))
        for bar, val in zip(bars, speedups):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                        f'{val:.1f}x', ha='center', va='bottom', fontsize=8)

    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='T1 baseline')
    ax.set_xlabel('Scenario')
    ax.set_ylabel('Speedup vs T1')
    ax.set_title('Speedup Relative to T1 Baseline')
    ax.set_xticks(list(x_base))
    ax.set_xticklabels(scenarios, rotation=15)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_dir / 'speedup.png', dpi=150)
    plt.close(fig)
    print(f'  -> {output_dir / "speedup.png"}')


def plot_throughput(results: dict, output_dir: Path):
    """Bar chart: throughput (reads/sec) per tier."""
    if 'throughput' not in results:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    tiers = []
    rates = []
    colors = []

    for tier in TIER_ORDER:
        data = results['throughput'].get(tier, {})
        if not data:
            continue
        median = data['median']
        n_reads = data.get('n_reads', 500)
        reads_per_sec = n_reads / median if median > 0 else 0
        tiers.append(TIER_LABELS.get(tier, tier))
        rates.append(reads_per_sec)
        colors.append(TIER_COLORS.get(tier, '#999'))

    bars = ax.bar(tiers, rates, color=colors)
    for bar, val in zip(bars, rates):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f'{val:,.0f}', ha='center', va='bottom', fontsize=9)

    ax.set_ylabel('Reads / second')
    ax.set_title('Throughput (500-read batch)')
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_dir / 'throughput.png', dpi=150)
    plt.close(fig)
    print(f'  -> {output_dir / "throughput.png"}')


def print_summary_table(results: dict):
    """Print a text summary table to stdout."""
    print('\n' + '=' * 72)
    print('BENCHMARK SUMMARY')
    print('=' * 72)

    if 'latency' in results:
        print('\nSingle-Read Latency (ms, median):')
        scenarios = list(results['latency'].keys())
        header = f'{"Scenario":<20s}'
        for tier in LATENCY_TIERS:
            header += f'{TIER_LABELS.get(tier, tier):>14s}'
        print(header)
        print('-' * len(header))
        for sc in scenarios:
            row = f'{sc:<20s}'
            t1_val = results['latency'].get(sc, {}).get(
                'T1_original', {}).get('median', 0)
            for ti, tier in enumerate(LATENCY_TIERS):
                val = results['latency'].get(sc, {}).get(
                    tier, {}).get('median', 0)
                ms = val * 1000
                if tier != 'T1_original' and t1_val > 0 and val > 0:
                    speedup = t1_val / val
                    row += f'{ms:>10.3f}ms ({speedup:.1f}x)'[:14]
                else:
                    row += f'{ms:>10.3f}ms'[:14]
                row = row.ljust(20 + 14 * (ti + 1))
            print(row)

    if 'throughput' in results:
        # Detect core count for fair batch comparison
        import multiprocessing
        n_cores = multiprocessing.cpu_count()

        t1_data = results['throughput'].get('T1_original', {})
        t1_rate = 0
        if t1_data and t1_data.get('median', 0) > 0:
            t1_n = t1_data.get('n_reads', 500)
            t1_rate = t1_n / t1_data['median']

        print(f'\nThroughput (reads/sec, {n_cores} logical cores):')
        for tier in TIER_ORDER:
            data = results['throughput'].get(tier, {})
            if data:
                n = data.get('n_reads', 500)
                med = data['median']
                rate = n / med if med > 0 else 0
                label = TIER_LABELS.get(tier, tier)
                if tier == 'T5_batch_rayon' and t1_rate > 0:
                    # Fair comparison: vs T1 × N cores
                    t1_parallel = t1_rate * n_cores
                    fair_speedup = rate / t1_parallel
                    print(
                        f'  {label:<16s}: '
                        f'{rate:>10,.0f} reads/sec  '
                        f'({rate / t1_rate:.0f}x vs T1×1, '
                        f'{fair_speedup:.1f}x vs T1×{n_cores})'
                    )
                elif t1_rate > 0 and tier != 'T1_original':
                    print(
                        f'  {label:<16s}: '
                        f'{rate:>10,.0f} reads/sec  '
                        f'({rate / t1_rate:.1f}x vs T1)'
                    )
                else:
                    print(
                        f'  {label:<16s}: '
                        f'{rate:>10,.0f} reads/sec'
                    )
    print()


def main():
    if not RESULTS_FILE.exists():
        print(f'ERROR: {RESULTS_FILE} not found.')
        print('Run benchmarks first:')
        print('  pipenv run python -m pytest benchmarks/ '
              '--benchmark-only --benchmark-json=benchmarks/data/results.json')
        sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data = load_results(RESULTS_FILE)
    results = parse_benchmarks(data)

    print('Generating benchmark charts...')
    plot_latency_bars(results, OUTPUT_DIR)
    plot_speedup(results, OUTPUT_DIR)
    plot_throughput(results, OUTPUT_DIR)
    print_summary_table(results)
    print('Done.')


if __name__ == '__main__':
    main()
