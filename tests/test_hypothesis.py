"""Hypothesis property-based tests for codon_align().

These tests verify structural invariants hold for randomly generated
sequences, complementing the exact-match tests in test_codon_alignment.py.
"""

from hypothesis import assume, given, settings
from hypothesis import strategies as st

from postalign.models.na_position import NAPosition
from postalign.processors.codon_alignment import (
    REFGAP,
    SEQGAP,
    codon_align,
)
from tests.conftest import make_sequence, seq_to_str

DEFAULT_GPS: dict[int, dict[tuple[int, int], int]] = {
    REFGAP: {},
    SEQGAP: {},
}

BASES = 'ACGT'


def _insert_gaps(seq: str, gap_positions: list[int], gap_size: int) -> str:
    """Insert gap_size dashes at each position in gap_positions."""
    result = list(seq)
    for offset, pos in enumerate(sorted(gap_positions)):
        idx = pos + offset * gap_size
        for _ in range(gap_size):
            result.insert(idx, '-')
    return ''.join(result)


# Strategy: generate a codon-aligned base sequence (length divisible by 3)
codon_seq = st.integers(min_value=2, max_value=30).flatmap(
    lambda n_codons: st.text(
        alphabet=BASES, min_size=n_codons * 3, max_size=n_codons * 3
    )
)


class TestInvariantsHypothesis:
    @given(base_seq=codon_seq)
    @settings(max_examples=200)
    def test_no_gap_passthrough(self, base_seq: str) -> None:
        """When neither ref nor seq has gaps, output == input."""
        ref = make_sequence(base_seq, header='ref', seqid=0)
        seq = make_sequence(base_seq, header='seq', seqid=0)
        max_pos = NAPosition.max_pos(ref.seqtext)
        ref_out, seq_out = codon_align(ref, seq, 30, 10, DEFAULT_GPS, 1, max_pos)
        assert seq_to_str(ref_out) == base_seq
        assert seq_to_str(seq_out) == base_seq

    @given(
        base_seq=codon_seq,
        gap_codon_idx=st.integers(min_value=0, max_value=100),
    )
    @settings(max_examples=200)
    def test_single_3bp_deletion_preserves_content(
        self,
        base_seq: str,
        gap_codon_idx: int,
    ) -> None:
        """A 3bp deletion: non-gap content is preserved in output."""
        n_codons = len(base_seq) // 3
        assume(n_codons >= 3)
        gap_codon_idx = gap_codon_idx % (n_codons - 1)  # avoid very end
        gap_pos = gap_codon_idx * 3

        # ref = full base sequence, seq = base with 3bp gap
        ref_str = base_seq
        seq_str = base_seq[:gap_pos] + '---' + base_seq[gap_pos + 3 :]

        ref = make_sequence(ref_str, header='ref', seqid=0)
        seq = make_sequence(seq_str, header='seq', seqid=0)
        max_pos = NAPosition.max_pos(ref.seqtext)
        ref_out, seq_out = codon_align(ref, seq, 30, 10, DEFAULT_GPS, 1, max_pos)

        r = seq_to_str(ref_out)
        s = seq_to_str(seq_out)
        assert len(r) == len(s)
        assert r.replace('-', '') == ref_str
        assert s.replace('-', '') == seq_str.replace('-', '')

    @given(
        base_seq=codon_seq,
        insert_codon_idx=st.integers(min_value=0, max_value=100),
        insert_bases=st.text(alphabet=BASES, min_size=3, max_size=3),
    )
    @settings(max_examples=200)
    def test_single_3bp_insertion_preserves_content(
        self,
        base_seq: str,
        insert_codon_idx: int,
        insert_bases: str,
    ) -> None:
        """A 3bp insertion: non-gap content is preserved in output."""
        n_codons = len(base_seq) // 3
        assume(n_codons >= 4)
        # Avoid boundary positions (0 and last) — the algorithm strips
        # boundary gaps, which changes the effective sequence length
        insert_codon_idx = 1 + (insert_codon_idx % (n_codons - 2))
        ins_pos = insert_codon_idx * 3

        # ref has gap where seq has insertion
        ref_str = base_seq[:ins_pos] + '---' + base_seq[ins_pos:]
        seq_str = base_seq[:ins_pos] + insert_bases + base_seq[ins_pos:]

        ref = make_sequence(ref_str, header='ref', seqid=0)
        seq = make_sequence(seq_str, header='seq', seqid=0)
        max_pos = NAPosition.max_pos(ref.seqtext)
        ref_out, seq_out = codon_align(ref, seq, 30, 10, DEFAULT_GPS, 1, max_pos)

        r = seq_to_str(ref_out)
        s = seq_to_str(seq_out)
        assert len(r) == len(s)
        assert r.replace('-', '') == ref_str.replace('-', '')
        assert s.replace('-', '') == seq_str

    @given(
        base_seq=codon_seq,
        gap_codon_idx=st.integers(min_value=0, max_value=100),
        window_size=st.integers(min_value=1, max_value=30),
    )
    @settings(max_examples=200)
    def test_window_size_preserves_content(
        self,
        base_seq: str,
        gap_codon_idx: int,
        window_size: int,
    ) -> None:
        """Varying window_size must still preserve non-gap content."""
        n_codons = len(base_seq) // 3
        assume(n_codons >= 3)
        gap_codon_idx = gap_codon_idx % (n_codons - 1)
        gap_pos = gap_codon_idx * 3

        ref_str = base_seq
        seq_str = base_seq[:gap_pos] + '---' + base_seq[gap_pos + 3 :]

        ref = make_sequence(ref_str, header='ref', seqid=0)
        seq = make_sequence(seq_str, header='seq', seqid=0)
        max_pos = NAPosition.max_pos(ref.seqtext)
        ref_out, seq_out = codon_align(
            ref, seq, 30, window_size, DEFAULT_GPS, 1, max_pos
        )

        r = seq_to_str(ref_out)
        s = seq_to_str(seq_out)
        assert len(r) == len(s)
        assert r.replace('-', '') == ref_str
        assert s.replace('-', '') == seq_str.replace('-', '')
