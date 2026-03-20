"""Benchmark fixtures — generate realistic codon alignment inputs."""
import random

import pytest

from postalign.models.na_position import NAPosition
from postalign.models.position_flag import PositionFlag
from postalign.models.sequence import Sequence
from postalign.processors.codon_alignment import REFGAP, SEQGAP


BASES = b'ACGT'
DEFAULT_GPS: dict[int, dict[tuple[int, int], int]] = {REFGAP: {}, SEQGAP: {}}


def _rand_seq(length: int, rng: random.Random) -> str:
    return ''.join(chr(BASES[rng.randint(0, 3)]) for _ in range(length))


def _insert_gaps(
    seq: str, n_gaps: int, gap_size: int, rng: random.Random
) -> str:
    """Insert *n_gaps* blocks of *gap_size* dashes at random positions."""
    s = list(seq)
    for _ in range(n_gaps):
        pos = rng.randint(0, len(s))
        s[pos:pos] = ['-'] * gap_size
    return ''.join(s)


def make_pair(
    length: int, n_ref_gaps: int, n_seq_gaps: int,
    gap_size: int, seed: int,
) -> tuple[Sequence, Sequence, dict]:
    """Build a ref/seq Sequence pair with controlled gaps."""
    rng = random.Random(seed)
    ref_bases = _rand_seq(length, rng)
    seq_bases = list(ref_bases)
    # mutate ~5% of bases
    for i in range(len(seq_bases)):
        if rng.random() < 0.05:
            seq_bases[i] = chr(BASES[rng.randint(0, 3)])
    seq_bases_str = ''.join(seq_bases)

    ref_str = _insert_gaps(ref_bases, n_ref_gaps, gap_size, rng)
    seq_str = _insert_gaps(seq_bases_str, n_seq_gaps, gap_size, rng)

    # pad to equal length
    diff = len(ref_str) - len(seq_str)
    if diff > 0:
        seq_str += '-' * diff
    elif diff < 0:
        ref_str += '-' * (-diff)

    def _make_seq(s: str, header: str) -> Sequence:
        nas = []
        pos = 1
        for ch in s:
            n = ord(ch)
            if ch == '-':
                nas.append(NAPosition(n, -1, PositionFlag.NONE))
            else:
                nas.append(NAPosition(n, pos, PositionFlag.NONE))
                pos += 1
        return Sequence(
            header=header, description='', seqtext=nas,
            seqid=0, seqtype=NAPosition, abs_seqstart=0,
            skip_invalid=True,
        )

    return _make_seq(ref_str, 'ref'), _make_seq(seq_str, 'seq'), DEFAULT_GPS


# ---- Parametrized input fixtures ----

SCENARIOS = {
    'short_simple':   dict(length=150,  n_ref_gaps=1, n_seq_gaps=1, gap_size=3),
    'medium_multi':   dict(length=300,  n_ref_gaps=2, n_seq_gaps=2, gap_size=3),
    'long_complex':   dict(length=1000, n_ref_gaps=3, n_seq_gaps=3, gap_size=6),
    'very_long':      dict(length=3000, n_ref_gaps=4, n_seq_gaps=4, gap_size=9),
}


@pytest.fixture(params=list(SCENARIOS.keys()))
def scenario(request):
    """Return (name, ref_seq, seq_seq, gps, ref_start, ref_end)."""
    name = request.param
    kw = SCENARIOS[name]
    ref, seq, gps = make_pair(**kw, seed=42)
    ref_start = 1
    ref_end = NAPosition.max_pos(ref.seqtext)
    return name, ref, seq, gps, ref_start, ref_end
