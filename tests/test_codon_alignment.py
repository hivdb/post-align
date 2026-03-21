"""API-level tests for codon_align().

Tests the public API only — no internal functions — so the test suite
remains valid across the Python and Rust backends.

The autouse fixture ``_codon_align_impl`` parametrizes every test in
this module over the pure-Python/Cython backend and the Rust-accelerated
backend, ensuring identical outputs.
"""

import sys

import pytest

from postalign.processors.codon_alignment import (
    REFGAP,
    SEQGAP,
    codon_align as _python_codon_align,
)
from postalign.processors.codon_alignment_rust import (
    codon_align as _rust_codon_align,
)

from tests.conftest import make_run_codon_align

# Default: Python (overridden per-test by the fixture below)
run_codon_align = make_run_codon_align(_python_codon_align)

_IMPLS = {
    'python': _python_codon_align,
    'rust': _rust_codon_align,
}


@pytest.fixture(params=['python', 'rust'], autouse=True)
def _codon_align_impl(request, monkeypatch):
    """Swap module-level run_codon_align to use the selected backend."""
    impl = _IMPLS[request.param]
    monkeypatch.setattr(
        sys.modules[__name__],
        'run_codon_align',
        make_run_codon_align(impl),
    )


# ---------------------------------------------------------------------------
# 1. No gaps — should return unchanged
# ---------------------------------------------------------------------------

class TestNoGaps:

    def test_identical_sequences(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGATGATG', 'ATGATGATGATG')
        assert ref_out == 'ATGATGATGATG'
        assert seq_out == 'ATGATGATGATG'

    def test_mismatched_bases_no_gaps(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGATGATG', 'ATGCCCATGATG')
        assert ref_out == 'ATGATGATGATG'
        assert seq_out == 'ATGCCCATGATG'


# ---------------------------------------------------------------------------
# 2. Single insertion (gap in ref = REFGAP)
# ---------------------------------------------------------------------------

class TestSingleInsertion:

    def test_3bp_insertion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATG---ATGATG', 'ATGATGCCAATGATG')
        assert ref_out == 'ATGATG---ATGATG'
        assert seq_out == 'ATGATGCCAATGATG'

    def test_6bp_insertion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATG------ATGATG', 'ATGATGCCACCAATGATG')
        assert ref_out == 'ATGATG------ATGATG'
        assert seq_out == 'ATGATGCCACCAATGATG'

    def test_1bp_insertion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATG-ATGATG', 'ATGATGCATGATG')
        assert ref_out == 'ATGATG-ATGATG'
        assert seq_out == 'ATGATGCATGATG'

    def test_2bp_insertion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATG--ATGATG', 'ATGATGCCATGATG')
        assert ref_out == 'ATGATG--ATGATG'
        assert seq_out == 'ATGATGCCATGATG'

    def test_4bp_insertion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATG----ATGATG', 'ATGATGCCCCATGATG')
        assert ref_out == 'ATGATG----ATGATG'
        assert seq_out == 'ATGATGCCCCATGATG'

    def test_5bp_insertion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATG-----ATGATG', 'ATGATGCCCCCATGATG')
        assert ref_out == 'ATGATG-----ATGATG'
        assert seq_out == 'ATGATGCCCCCATGATG'


# ---------------------------------------------------------------------------
# 3. Single deletion (gap in seq = SEQGAP)
# ---------------------------------------------------------------------------

class TestSingleDeletion:

    def test_3bp_deletion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGCCAATGATG', 'ATGATG---ATGATG')
        assert ref_out == 'ATGATGCCAATGATG'
        assert seq_out == 'ATGATG---ATGATG'

    def test_6bp_deletion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGCCACCAATGATG', 'ATGATG------ATGATG')
        assert ref_out == 'ATGATGCCACCAATGATG'
        assert seq_out == 'ATGATG------ATGATG'

    def test_1bp_deletion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGCCAATGATG', 'ATGATG-CAATGATG')
        assert ref_out == 'ATGATGCCAATGATG'
        assert seq_out == 'ATGATGCA-ATGATG'

    def test_2bp_deletion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGCCAATGATG', 'ATGATG--AATGATG')
        assert ref_out == 'ATGATGCCAATGATG'
        assert seq_out == 'ATGATGA--ATGATG'

    def test_4bp_deletion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGCCAATGATGATG', 'ATGATG----ATGATGATG')
        assert ref_out == 'ATGATGCCAATGATGATG'
        assert seq_out == 'ATGATGATGATG---AT-G'

    def test_5bp_deletion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGCCAATGATGATGATG', 'ATGATG-----ATGATGATGATG')
        assert ref_out == 'ATGATGCCAATGATGATGATG'
        assert seq_out == 'ATGATGATGATGATG---A--TG'


# ---------------------------------------------------------------------------
# 4. Scattered gaps — core codon alignment behavior
#    Small gaps near each other should be gathered into one larger gap
# ---------------------------------------------------------------------------

class TestScatteredGaps:

    def test_two_1bp_deletions_gathered(self) -> None:
        """Two 1bp gaps near each other → gathered into one 2bp gap."""
        ref_out, seq_out = run_codon_align(
            'ATGATGCCAATGATGATG', 'ATGATG-C-ATGATGATG')
        assert ref_out == 'ATGATGCCAATGATGATG'
        assert seq_out == 'ATGATGC--ATGATGATG'

    def test_three_1bp_deletions_gathered(self) -> None:
        """Three 1bp gaps → gathered into one 3bp gap."""
        ref_out, seq_out = run_codon_align(
            'ATGATGCCAATGATGATG', 'ATGATG-C-A-GATGATG')
        assert ref_out == 'ATGATGCCAATGATGATG'
        assert seq_out == 'ATGATGCAGATGATG---'

    def test_2bp_plus_1bp_deletions_gathered(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGCCAATGATGATG', 'ATGATG--A-TGATGATG')
        assert ref_out == 'ATGATGCCAATGATGATG'
        assert seq_out == 'ATGATG---ATGATGATG'

    def test_1bp_plus_2bp_deletions_gathered(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGCCAATGATGATG', 'ATGATG-CA--GATGATG')
        assert ref_out == 'ATGATGCCAATGATGATG'
        assert seq_out == 'ATGATGCAGATGATG---'

    def test_two_2bp_deletions_gathered(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGCCAAATGATGATG', 'ATGATG--A--GATGATG')
        assert ref_out == 'ATGATGCCAAATGATGATG'
        assert seq_out == 'ATGATGAGATGA---TG-'

    def test_two_1bp_insertions_gathered(self) -> None:
        """Two 1bp ref gaps near each other → gathered."""
        ref_out, seq_out = run_codon_align(
            'ATGATG-C-ATGATG', 'ATGATGCCAATGATG')
        assert ref_out == 'ATGATG--CATGATG'
        assert seq_out == 'ATGATGCCAATGATG'

    def test_three_1bp_insertions_gathered(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATG-C-A-GATG', 'ATGATGCCAAGATG')
        assert ref_out == 'ATGATG---CAGATG'
        assert seq_out == 'ATGATGCCAAGATG'

    def test_2bp_plus_1bp_insertions_gathered(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATG--A-TGATG', 'ATGATGCCAATGATG')
        assert ref_out == 'ATGATG---ATGATG'
        assert seq_out == 'ATGATGCCAATGATG'

    def test_1bp_plus_2bp_insertions_gathered(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATG-CA--GATG', 'ATGATGCCAAGATG')
        assert ref_out == 'ATGATG---CAGATG'
        assert seq_out == 'ATGATGCCAAGATG'


# ---------------------------------------------------------------------------
# 5. Multiple gaps (triplet)
# ---------------------------------------------------------------------------

class TestMultipleGaps:

    def test_insertion_and_deletion(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATG---CCAATG---ATG',
            'ATGATGCCACCA---ATGATG')
        assert ref_out == 'ATGATG---CCAATGATG'
        assert seq_out == 'ATGATGCCACCAATGATG'

    def test_multiple_insertions(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATG---ATG---ATG', 'ATGCCAATGGGGATG')
        assert ref_out == 'ATG------ATGATG'
        assert seq_out == 'ATGCCAATGGGGATG'

    def test_multiple_deletions(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGCCAATGGGGATG', 'ATG---ATG---ATG')
        assert ref_out == 'ATGCCAATGGGGATG'
        assert seq_out == '------ATGATGATG'


# ---------------------------------------------------------------------------
# 6. Boundary gaps
# ---------------------------------------------------------------------------

class TestBoundaryGaps:

    def test_gap_at_sequence_start(self) -> None:
        ref_out, seq_out = run_codon_align(
            '---ATGATGATG', 'CCAATGATGATG')
        assert ref_out == 'ATGATGATG'
        assert seq_out == 'ATGATGATG'

    def test_gap_at_sequence_end(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGATG---', 'ATGATGATGCCA')
        assert ref_out == 'ATGATGATG---'
        assert seq_out == 'ATGATGATGCCA'

    def test_deletion_at_start(self) -> None:
        ref_out, seq_out = run_codon_align(
            'CCAATGATGATG', '---ATGATGATG')
        assert ref_out == 'CCAATGATGATG'
        assert seq_out == '---ATGATGATG'

    def test_deletion_at_end(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGATGCCA', 'ATGATGATG---')
        assert ref_out == 'ATGATGATGCCA'
        assert seq_out == 'ATGATGATG---'


# ---------------------------------------------------------------------------
# 7. Adjacent/separated gaps (min_gap_distance)
# ---------------------------------------------------------------------------

class TestGapDistance:

    def test_adjacent_gaps_within_distance(self) -> None:
        """Gaps within min_gap_distance should be gathered into one window."""
        ref_out, seq_out = run_codon_align(
            'ATGATGATGATGATGATGATGATGATGATG',
            'ATG---ATGATGATGATG---ATGATGATG',
            min_gap_distance=30)
        assert ref_out == 'ATGATGATGATGATGATGATGATGATGATG'
        assert seq_out == '------ATGATGATGATGATGATGATGATG'

    def test_separated_gaps_beyond_distance(self) -> None:
        """Gaps beyond min_gap_distance should be in separate windows."""
        ref_out, seq_out = run_codon_align(
            'ATGATGATGATGATGATGATGATGATGATGATGATGATGATG',
            'ATG---ATGATGATGATGATGATGATGATGATGATG---ATG',
            min_gap_distance=6)
        assert ref_out == 'ATGATGATGATGATGATGATGATGATGATGATGATGATGATG'
        assert seq_out == '---ATGATGATGATGATGATGATGATGATGATGATGATG---'


# ---------------------------------------------------------------------------
# 8. Ambiguous bases (IUPAC)
# ---------------------------------------------------------------------------

class TestAmbiguousBases:

    @pytest.mark.parametrize('ambig', [
        'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N',
    ])
    def test_iupac_codes_preserve_content(self, ambig: str) -> None:
        ref = f'ATG{ambig}{ambig}{ambig}ATGATG'
        seq = 'ATG---ATGATG'
        ref_out, seq_out = run_codon_align(ref, seq)
        assert len(ref_out) == len(seq_out) == len(ref)
        assert ref_out.replace('-', '') == ref.replace('-', '')
        assert seq_out.replace('-', '') == seq.replace('-', '')


# ---------------------------------------------------------------------------
# 9. Short / edge-case sequences
# ---------------------------------------------------------------------------

class TestEdgeCases:

    def test_short_sequence_less_than_codon(self) -> None:
        ref_out, seq_out = run_codon_align('AT', 'AT')
        assert ref_out == 'AT'
        assert seq_out == 'AT'

    def test_single_codon(self) -> None:
        ref_out, seq_out = run_codon_align('ATG', 'ATG')
        assert ref_out == 'ATG'
        assert seq_out == 'ATG'

    def test_all_gaps_in_seq(self) -> None:
        ref_out, seq_out = run_codon_align('ATGATGATG', '---------')
        assert ref_out == 'ATGATGATG'
        assert seq_out == '---------'


# ---------------------------------------------------------------------------
# 10. Redundant gap removal
# ---------------------------------------------------------------------------

class TestRedundantGaps:

    def test_matched_gaps_removed(self) -> None:
        """When both ref and seq have gaps at same region, redundant
        gaps should be removed (sequence gets shorter)."""
        ref_out, seq_out = run_codon_align(
            'ATG---ATGATG', 'ATG---ATGATG')
        assert ref_out == 'ATGATGATG'
        assert seq_out == 'ATGATGATG'


# ---------------------------------------------------------------------------
# 11. ref_start / ref_end (partial range)
# ---------------------------------------------------------------------------

class TestPartialRange:

    def test_ref_start_limits_application(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATG---ATGATGATGATG', 'ATGCCAATGATGATGATG',
            ref_start=4)
        assert ref_out == 'ATGATGATGATGATG'
        assert seq_out == 'ATGATGATGATGATG'

    def test_ref_end_limits_application(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGATG---ATGATG', 'ATGATGATGCCAATGATG',
            ref_end=6)
        assert ref_out == 'ATGATGATG---ATGATG'
        assert seq_out == 'ATGATGATGCCAATGATG'


# ---------------------------------------------------------------------------
# 12. Custom gap_placement_score
# ---------------------------------------------------------------------------

class TestGapPlacementScore:

    def test_bonus_score_keeps_gap_in_place(self) -> None:
        """A large bonus at position 7 keeps the gap at position 7."""
        gps: dict[int, dict[tuple[int, int], int]] = {
            REFGAP: {},
            SEQGAP: {(7, 0): 100},
        }
        ref_out, seq_out = run_codon_align(
            'ATGATGATGATGATGATG', 'ATGATG---ATGATGATG',
            gap_placement_score=gps)
        assert ref_out == 'ATGATGATGATGATGATG'
        assert seq_out == 'ATGATG---ATGATGATG'

    def test_penalty_score_moves_gap(self) -> None:
        """A large penalty at position 7 pushes the gap elsewhere."""
        gps: dict[int, dict[tuple[int, int], int]] = {
            REFGAP: {},
            SEQGAP: {(7, 0): -100},
        }
        ref_out, seq_out = run_codon_align(
            'ATGATGATGATGATGATG', 'ATGATG---ATGATGATG',
            gap_placement_score=gps)
        assert ref_out == 'ATGATGATGATGATGATG'
        assert seq_out == '---ATGATGATGATGATG'


# ---------------------------------------------------------------------------
# 13. Window size parameter
# ---------------------------------------------------------------------------

class TestWindowSize:

    def test_small_window(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGATGATGATGATG', 'ATGATG---ATGATGATG',
            window_size=1)
        assert ref_out == 'ATGATGATGATGATGATG'
        assert seq_out == 'ATG---ATGATGATGATG'

    def test_large_window(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGATGATGATGATGATG', 'ATGATG---ATGATGATG',
            window_size=50)
        assert ref_out == 'ATGATGATGATGATGATG'
        assert seq_out == '---ATGATGATGATGATG'


# ---------------------------------------------------------------------------
# 14. Invariants (applicable to all cases)
# ---------------------------------------------------------------------------

class TestInvariants:
    """Cross-cutting invariants that must hold for any codon_align call."""

    @pytest.mark.parametrize('ref,seq', [
        ('ATGATGATGATG', 'ATGATGATGATG'),
        ('ATGATG---ATGATG', 'ATGATGCCAATGATG'),
        ('ATGATGCCAATGATG', 'ATGATG---ATGATG'),
        ('ATG---ATG---ATG', 'ATGCCAATGGGGATG'),
        ('---ATGATGATG', 'CCAATGATGATG'),
        ('ATGATGATG---', 'ATGATGATGCCA'),
    ])
    def test_output_lengths_match(self, ref: str, seq: str) -> None:
        ref_out, seq_out = run_codon_align(ref, seq)
        assert len(ref_out) == len(seq_out)

    @pytest.mark.parametrize('ref,seq', [
        ('ATGATG---ATGATG', 'ATGATGCCAATGATG'),
        ('ATGATGCCAATGATG', 'ATGATG---ATGATG'),
        ('ATG---ATG---ATG', 'ATGCCAATGGGGATG'),
    ])
    def test_nongap_content_preserved(self, ref: str, seq: str) -> None:
        ref_out, seq_out = run_codon_align(ref, seq)
        assert ref_out.replace('-', '') == ref.replace('-', '')
        assert seq_out.replace('-', '') == seq.replace('-', '')


# ===================================================================
# Edge case tests (A–I) — realistic NA-NA aligner artifacts
# ===================================================================

# Heterogeneous codon sequences for tests where BLOSUM62 matters
# HETERO = Met-Pro-Glu-Trp-Cys-Lys (6 codons, 18bp)
HETERO = 'ATGCCTGAATGCTGTAAG'
# HETERO9 = 9 codons, 27bp
HETERO9 = 'ATGCCTGAATGCTGTAAGTTTGCATAT'
# LONG_REF = 12 codons, 36bp
LONG_REF = 'ATGCCTGAATGCTGTAAGTTTGCATATCAGACTGTT'


# ---------------------------------------------------------------------------
# A. Frameshift correction — gaps placed off codon boundary by NA aligner
#    Uses heterogeneous codons so BLOSUM62 can differentiate positions.
# ---------------------------------------------------------------------------

class TestFrameshiftCorrection:

    def test_3bp_del_shifted_plus1(self) -> None:
        """3bp deletion gap shifted +1 into CCT codon."""
        ref_out, seq_out = run_codon_align(
            HETERO, 'ATGC---AATGCTGTAAG')
        assert ref_out == 'ATGCCTGAATGCTGTAAG'
        assert seq_out == 'ATG---CAATGCTGTAAG'

    def test_3bp_del_shifted_plus2(self) -> None:
        """3bp deletion gap shifted +2 into CCT codon."""
        ref_out, seq_out = run_codon_align(
            HETERO, 'ATGCC---ATGCTGTAAG')
        assert ref_out == 'ATGCCTGAATGCTGTAAG'
        assert seq_out == 'ATGCCA---TGCTGTAAG'

    def test_6bp_del_shifted_plus1(self) -> None:
        """6bp deletion shifted +1 off codon boundary."""
        ref_out, seq_out = run_codon_align(
            HETERO9, 'ATGC------TGCTGTAAGTTTGCATAT')
        assert ref_out == 'ATGCCTGAATGCTGTAAGTTTGCATAT'
        assert seq_out == 'ATG------CTGCTGTAAGTTTGCATAT'

    def test_6bp_del_shifted_plus2(self) -> None:
        """6bp deletion shifted +2 off codon boundary."""
        ref_out, seq_out = run_codon_align(
            HETERO9, 'ATGCC------GCTGTAAGTTTGCATAT')
        assert ref_out == 'ATGCCTGAATGCTGTAAGTTTGCATAT'
        assert seq_out == 'ATGCCG------CTGTAAGTTTGCATAT'

    def test_3bp_ins_shifted_plus1(self) -> None:
        """3bp insertion with ref gap shifted +1 off boundary."""
        ref_out, seq_out = run_codon_align(
            'ATGC---CTGAATGCTGTAAG',
            'ATGCAAACTGAATGCTGTAAG')
        assert ref_out == 'ATG---CCTGAATGCTGTAAG'
        assert seq_out == 'ATGCAAACTGAATGCTGTAAG'

    def test_1bp_del_frameshift_codon_pos1(self) -> None:
        """1bp deletion (true frameshift) splitting codon at pos 1.
        Gap ends up at codon end due to move_gap_to_codon_end."""
        ref_out, seq_out = run_codon_align(
            HETERO, 'ATGC-TGAATGCTGTAAG')
        assert ref_out == 'ATGCCTGAATGCTGTAAG'
        assert seq_out == 'ATGCT-GAATGCTGTAAG'

    def test_1bp_del_frameshift_codon_pos2(self) -> None:
        """1bp deletion (true frameshift) splitting codon at pos 2."""
        ref_out, seq_out = run_codon_align(
            HETERO, 'ATGCC-GAATGCTGTAAG')
        assert ref_out == 'ATGCCTGAATGCTGTAAG'
        assert seq_out == 'ATGCCGAA-TGCTGTAAG'

    def test_2bp_del_frameshift(self) -> None:
        """2bp deletion (frameshift) off boundary."""
        ref_out, seq_out = run_codon_align(
            HETERO, 'ATGC--GAATGCTGTAAG')
        assert ref_out == 'ATGCCTGAATGCTGTAAG'
        assert seq_out == 'ATGC--GAATGCTGTAAG'


# ---------------------------------------------------------------------------
# B. Maximally scattered gaps — worst-case noisy aligner output
# ---------------------------------------------------------------------------

class TestMaximallyScatteredGaps:

    def test_every_other_base_deletion(self) -> None:
        """Every other base is a gap — all should gather."""
        ref_out, seq_out = run_codon_align(
            'ATGCCTGAATGCTGT',
            'A-G-C-G-A-G-T-T')
        assert ref_out == 'ATGCCTGAATGCTGT'
        assert seq_out == '------AG-CGAGTT'

    def test_scattered_partial(self) -> None:
        """Scattered gaps in first half, normal second half."""
        ref_out, seq_out = run_codon_align(
            HETERO,
            'A-G-C-G-A-GCTGTAAG')
        assert ref_out == 'ATGCCTGAATGCTGTAAG'
        assert seq_out == '---A--GCGAGCTGTAAG'

    def test_scattered_insertions_in_ref(self) -> None:
        """Scattered 1bp ref gaps (insertions) gathered."""
        ref_out, seq_out = run_codon_align(
            'ATG-C-C-TGAATGC',
            'ATGACCCCTGAATGC')
        assert ref_out == 'ATG---CCTGAATGC'
        assert seq_out == 'ATGACCCCTGAATGC'


# ---------------------------------------------------------------------------
# C. Homopolymer gaps — aligners struggle with gap placement in runs
# ---------------------------------------------------------------------------

class TestHomopolymerGaps:

    def test_polyA_deletion(self) -> None:
        """3bp deletion in poly-A run, already on boundary."""
        ref_out, seq_out = run_codon_align(
            'ATGAAAAAGATGATG',
            'ATGAA---GATGATG')
        assert ref_out == 'ATGAAAAAGATGATG'
        assert seq_out == 'ATG---AAGATGATG'

    def test_polyT_deletion(self) -> None:
        """3bp deletion in poly-T run."""
        ref_out, seq_out = run_codon_align(
            'ATGTTTTTTATGATG',
            'ATGTTT---ATGATG')
        assert ref_out == 'ATGTTTTTTATGATG'
        assert seq_out == 'ATG---TTTATGATG'

    def test_polyA_off_boundary(self) -> None:
        """3bp deletion in poly-A, gap placed off boundary by aligner."""
        ref_out, seq_out = run_codon_align(
            'ATGAAAAAGATGATG',
            'ATGA---AGATGATG')
        assert ref_out == 'ATGAAAAAGATGATG'
        assert seq_out == 'ATG---AAGATGATG'

    def test_scattered_in_homopolymer(self) -> None:
        """Three 1bp gaps scattered in poly-A → gathered."""
        ref_out, seq_out = run_codon_align(
            'ATGAAAAAGATGATG',
            'ATG-A-A-GATGATG')
        assert ref_out == 'ATGAAAAAGATGATG'
        assert seq_out == 'ATG---AAGATGATG'


# ---------------------------------------------------------------------------
# D. Gap adjacent to mismatch (SNP) — common aligner artifact
# ---------------------------------------------------------------------------

class TestGapAdjacentToMismatch:

    def test_3bp_del_next_to_snp(self) -> None:
        """3bp deletion with a SNP in the adjacent codon."""
        ref_out, seq_out = run_codon_align(
            'ATGCCTGAAATGATG',
            'ATGCCT---CTGATG')
        assert ref_out == 'ATGCCTGAAATGATG'
        assert seq_out == 'ATGCCT---CTGATG'

    def test_3bp_del_with_multiple_snps(self) -> None:
        """3bp deletion with SNPs in multiple codons."""
        ref_out, seq_out = run_codon_align(
            'ATGCCTGAAATGATG',
            'CTGCCT---CTGATG')
        assert ref_out == 'ATGCCTGAAATGATG'
        assert seq_out == 'CTGCCT---CTGATG'

    def test_1bp_del_next_to_snp(self) -> None:
        """1bp frameshift deletion adjacent to a SNP."""
        ref_out, seq_out = run_codon_align(
            'ATGCCTGAAATGATG',
            'ATGCCT-GAATGATG')
        assert ref_out == 'ATGCCTGAAATGATG'
        assert seq_out == 'ATGCCTGAATG-ATG'


# ---------------------------------------------------------------------------
# E. Long deletions (multi-codon) — HIV drug resistance regions
# ---------------------------------------------------------------------------

class TestLongDeletions:

    def test_9bp_deletion_on_boundary(self) -> None:
        """9bp (3-codon) deletion already on boundary."""
        ref_out, seq_out = run_codon_align(
            LONG_REF,
            'ATG---------TGTAAGTTTGCATATCAGACTGTT')
        assert ref_out == 'ATGCCTGAATGCTGTAAGTTTGCATATCAGACTGTT'
        assert seq_out == 'ATG---------TGTAAGTTTGCATATCAGACTGTT'

    def test_12bp_deletion_on_boundary(self) -> None:
        """12bp (4-codon) deletion."""
        ref_out, seq_out = run_codon_align(
            LONG_REF,
            'ATG------------AAGTTTGCATATCAGACTGTT')
        assert ref_out == 'ATGCCTGAATGCTGTAAGTTTGCATATCAGACTGTT'
        assert seq_out == '------------ATGAAGTTTGCATATCAGACTGTT'

    def test_9bp_deletion_off_boundary(self) -> None:
        """9bp deletion shifted +1 off codon boundary."""
        ref_out, seq_out = run_codon_align(
            LONG_REF,
            'ATGC---------GTAAGTTTGCATATCAGACTGTT')
        assert ref_out == 'ATGCCTGAATGCTGTAAGTTTGCATATCAGACTGTT'
        assert seq_out == 'ATGCGT---------AAGTTTGCATATCAGACTGTT'


# ---------------------------------------------------------------------------
# F. Mixed insertion + deletion in same window
# ---------------------------------------------------------------------------

class TestMixedInsDelSameWindow:

    def test_ref_and_seq_gap_same_window(self) -> None:
        """Both ref and seq have 3bp gap in same window → redundant
        gaps removed, remaining placed optimally."""
        ref_out, seq_out = run_codon_align(
            'ATGCCT---GAATGCTGT',
            'ATGCCTGAA---TGCTGT')
        assert ref_out == 'ATGCCTGAATGCTGT'
        assert seq_out == 'ATGCCTGAATGCTGT'

    def test_overlapping_ins_del(self) -> None:
        """Ref has insertion gap, seq has deletion gap at different
        positions → redundant removal + alignment."""
        ref_out, seq_out = run_codon_align(
            'ATG---CCTGAATGCTGT',
            'ATGCCTGAA---TGCTGT')
        assert ref_out == 'ATGCCTGAATGCTGT'
        assert seq_out == 'ATGCCTGAATGCTGT'

    def test_unequal_ins_del(self) -> None:
        """Ref has 3bp gap, seq has 6bp gap → net 3bp deletion."""
        ref_out, seq_out = run_codon_align(
            'ATG---CCTGAATGCTGTAAG',
            'ATGCCTGAATGC------AAG')
        assert ref_out == 'ATGCCTGAATGCTGTAAG'
        assert seq_out == 'ATGCCTGAATGC---AAG'


# ---------------------------------------------------------------------------
# G. Gap placement score interacting with frameshift
# ---------------------------------------------------------------------------

class TestGapPlacementScoreWithFrameshift:

    def test_bonus_at_pos4(self) -> None:
        """Bonus score at position 4 influences gap placement."""
        gps: dict[int, dict[tuple[int, int], int]] = {
            REFGAP: {},
            SEQGAP: {(4, 0): 100},
        }
        ref_out, seq_out = run_codon_align(
            HETERO, 'ATGCC---ATGCTGTAAG',
            gap_placement_score=gps)
        assert ref_out == 'ATGCCTGAATGCTGTAAG'
        assert seq_out == 'ATG---CCATGCTGTAAG'

    def test_penalty_at_pos4(self) -> None:
        """Penalty at position 4 pushes gap away."""
        gps: dict[int, dict[tuple[int, int], int]] = {
            REFGAP: {},
            SEQGAP: {(4, 0): -100},
        }
        ref_out, seq_out = run_codon_align(
            HETERO, 'ATGCC---ATGCTGTAAG',
            gap_placement_score=gps)
        assert ref_out == 'ATGCCTGAATGCTGTAAG'
        assert seq_out == 'ATGCCA---TGCTGTAAG'

    def test_size_specific_score(self) -> None:
        """Size-specific score (pos, 3) takes precedence over (pos, 0)."""
        gps: dict[int, dict[tuple[int, int], int]] = {
            REFGAP: {},
            SEQGAP: {(4, 3): 100, (4, 0): -100},
        }
        ref_out, seq_out = run_codon_align(
            HETERO, 'ATGCCT---GCTGTAAG',
            gap_placement_score=gps)
        assert ref_out == 'ATGCCTGAATGCTGTAAG'
        assert seq_out == 'ATG---CCTGCTGTAAG'

    def test_competing_scores(self) -> None:
        """Competing bonus/penalty at different positions."""
        gps: dict[int, dict[tuple[int, int], int]] = {
            REFGAP: {},
            SEQGAP: {(4, 0): 100, (10, 0): -100},
        }
        ref_out, seq_out = run_codon_align(
            HETERO, 'ATGCCT---GCTGTAAG',
            gap_placement_score=gps)
        assert ref_out == 'ATGCCTGAATGCTGTAAG'
        assert seq_out == 'ATG---CCTGCTGTAAG'


# ---------------------------------------------------------------------------
# H. Adversarial / pathological inputs
# ---------------------------------------------------------------------------

class TestAdversarialInputs:

    def test_seq_entirely_gaps(self) -> None:
        ref_out, seq_out = run_codon_align('ATGATGATG', '---------')
        assert ref_out == 'ATGATGATG'
        assert seq_out == '---------'

    def test_both_entirely_gaps(self) -> None:
        ref_out, seq_out = run_codon_align('---------', '---------')
        assert ref_out == '---------'
        assert seq_out == '---------'

    def test_single_base_among_gaps(self) -> None:
        ref_out, seq_out = run_codon_align('ATGATGATG', '----A----')
        assert ref_out == 'ATGATGATG'
        assert seq_out == '------A--'

    def test_gap_longer_than_content(self) -> None:
        ref_out, seq_out = run_codon_align(
            'ATGCCTGAATGCTGT', 'A--------------')
        assert ref_out == 'ATGCCTGAATGCTGT'
        assert seq_out == 'A--------------'

    def test_long_sequence_150bp(self) -> None:
        """100+ codons with a single 3bp gap — no crash."""
        ref = 'ATG' * 50
        seq = 'ATG' * 16 + '---' + 'ATG' * 33
        ref_out, seq_out = run_codon_align(ref, seq)
        assert len(ref_out) == len(seq_out)
        assert ref_out.replace('-', '') == ref
        assert seq_out.replace('-', '') == seq.replace('-', '')

    def test_unequal_gap_counts(self) -> None:
        """Ref has 3 gaps, seq has 6 gaps → net deletion after
        redundant removal."""
        ref_out, seq_out = run_codon_align(
            'ATG---CCTGAATGC',
            'ATGCCT------TGC')
        assert ref_out == 'ATGCCTGAATGC'
        assert seq_out == 'ATGCCT---TGC'


# ---------------------------------------------------------------------------
# I. Redundant gap edge cases
# ---------------------------------------------------------------------------

class TestRedundantGapEdgeCases:

    def test_asymmetric_ref3_seq1(self) -> None:
        """Ref has 3 gaps, seq has 1 gap → 1 removed from each,
        net 2 ref gaps remain."""
        ref_out, seq_out = run_codon_align(
            'ATG---CCTGAATGC',
            'ATG-CCTGAATGCCC')
        assert ref_out == 'ATGCCTGAATGC--'
        assert seq_out == 'ATGCCTGAATGCCC'

    def test_redundant_at_different_positions(self) -> None:
        """Ref and seq both have 3bp gaps at different positions →
        all redundant, both shrink."""
        ref_out, seq_out = run_codon_align(
            'ATG---CCTGAATGCTGT',
            'ATGCCT---GAATGCTGT')
        assert ref_out == 'ATGCCTGAATGCTGT'
        assert seq_out == 'ATGCCTGAATGCTGT'

    def test_large_redundant_6bp(self) -> None:
        """Both have 6bp gaps at same position → all removed."""
        ref_out, seq_out = run_codon_align(
            'ATG------GAATGCTGT',
            'ATG------GAATGCTGT')
        assert ref_out == 'ATGGAATGCTGT'
        assert seq_out == 'ATGGAATGCTGT'
