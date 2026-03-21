"""Tests for internal/helper functions in codon_alignment.py.

Targets functions not exercised through the public codon_align() API
to boost coverage toward >95%.
"""

import click
import pytest

from postalign.models.na_position import NAPosition
from postalign.models.position_flag import PositionFlag
from postalign.processors.codon_alignment import (
    LEFT,
    NOGAP,
    REFGAP,
    RIGHT,
    SEQGAP,
    calc_match_score_precomputed,
    center_expand_positions,
    extend_codons_until_gap,
    find_first_gap,
    find_windows_with_gap,
    gap_placement_score_callback,
    gather_gaps,
    move_gap_to_codon_end,
    parse_gap_placement_score,
    realign_gaps,
    remove_n_gaps,
    remove_redundant_gaps,
    separate_gaps_from_nas,
)
from postalign.utils.codonutils import translate_codons
from tests.conftest import make_na_positions

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _na(ch: str, pos: int = -1) -> NAPosition:
    """Shorthand for creating a single NAPosition."""
    return NAPosition(ord(ch), pos, PositionFlag.NONE)


def _nas(s: str) -> list[NAPosition]:
    """Create NAPositions from a string (gaps get pos=-1)."""
    return make_na_positions(s)


# ---------------------------------------------------------------------------
# find_first_gap
# ---------------------------------------------------------------------------


class TestFindFirstGap:
    def test_no_gap(self) -> None:
        assert find_first_gap(_nas('ATGATG')) == -1

    def test_gap_at_start(self) -> None:
        assert find_first_gap(_nas('-ATGATG')) == 0

    def test_gap_in_middle(self) -> None:
        assert find_first_gap(_nas('ATG---ATG')) == 3

    def test_gap_at_end(self) -> None:
        assert find_first_gap(_nas('ATG-')) == 3

    def test_empty_list(self) -> None:
        assert find_first_gap([]) == -1

    def test_dot_gap(self) -> None:
        assert find_first_gap(_nas('ATG.ATG')) == 3


# ---------------------------------------------------------------------------
# separate_gaps_from_nas
# ---------------------------------------------------------------------------


class TestSeparateGapsFromNas:
    def test_no_gaps(self) -> None:
        nas = _nas('ATGATG')
        nongaps, gaps = separate_gaps_from_nas(nas)
        assert len(nongaps) == 6
        assert len(gaps) == 0

    def test_all_gaps(self) -> None:
        nas = _nas('---')
        nongaps, gaps = separate_gaps_from_nas(nas)
        assert len(nongaps) == 0
        assert len(gaps) == 3

    def test_mixed(self) -> None:
        nas = _nas('A-G-T')
        nongaps, gaps = separate_gaps_from_nas(nas)
        assert len(nongaps) == 3
        assert len(gaps) == 2
        assert all(not n.is_gap for n in nongaps)
        assert all(n.is_gap for n in gaps)


# ---------------------------------------------------------------------------
# remove_n_gaps
# ---------------------------------------------------------------------------


class TestRemoveNGaps:
    def test_remove_zero(self) -> None:
        nas = _nas('A-G-T')
        result = remove_n_gaps(nas, 0)
        assert len(result) == 5

    def test_remove_one(self) -> None:
        nas = _nas('A-G-T')
        result = remove_n_gaps(nas, 1)
        assert len(result) == 4

    def test_remove_all_gaps(self) -> None:
        nas = _nas('A-G-T')
        result = remove_n_gaps(nas, 2)
        assert len(result) == 3
        assert all(not n.is_gap for n in result)

    def test_remove_more_than_available(self) -> None:
        nas = _nas('A-T')
        result = remove_n_gaps(nas, 5)
        assert len(result) == 2
        assert all(not n.is_gap for n in result)


# ---------------------------------------------------------------------------
# remove_redundant_gaps
# ---------------------------------------------------------------------------


class TestRemoveRedundantGaps:
    def test_no_redundant(self) -> None:
        ref = _nas('ATG---ATG')
        seq = _nas('ATGCCCATG')
        r, s = remove_redundant_gaps(ref, seq)
        assert len(r) == 9
        assert len(s) == 9

    def test_equal_gaps(self) -> None:
        ref = _nas('ATG---ATG')
        seq = _nas('ATG---ATG')
        r, s = remove_redundant_gaps(ref, seq)
        assert NAPosition.count_gaps(r) == 0
        assert NAPosition.count_gaps(s) == 0

    def test_unequal_gaps(self) -> None:
        ref = _nas('A---TG')
        seq = _nas('A-CCTG')
        r, s = remove_redundant_gaps(ref, seq)
        ref_gaps = NAPosition.count_gaps(r)
        seq_gaps = NAPosition.count_gaps(s)
        assert min(ref_gaps, seq_gaps) == 0


# ---------------------------------------------------------------------------
# move_gap_to_codon_end
# ---------------------------------------------------------------------------


class TestMoveGapToCodonEnd:
    def test_gap_already_at_end(self) -> None:
        codon = [_na('A', 1), _na('T', 2), _na('-')]
        result = move_gap_to_codon_end([codon])
        assert result[0][-1].is_gap

    def test_gap_at_start(self) -> None:
        codon = [_na('-'), _na('A', 1), _na('T', 2)]
        result = move_gap_to_codon_end([codon])
        assert not result[0][0].is_gap
        assert result[0][-1].is_gap

    def test_no_gaps(self) -> None:
        codon = [_na('A', 1), _na('T', 2), _na('G', 3)]
        result = move_gap_to_codon_end([codon])
        assert len(result[0]) == 3
        assert not any(n.is_gap for n in result[0])

    def test_multiple_codons(self) -> None:
        c1 = [_na('-'), _na('A', 1), _na('T', 2)]
        c2 = [_na('G', 3), _na('-'), _na('C', 4)]
        result = move_gap_to_codon_end([c1, c2])
        for codon in result:
            nongaps = [n for n in codon if not n.is_gap]
            gaps = [n for n in codon if n.is_gap]
            assert codon == nongaps + gaps


# ---------------------------------------------------------------------------
# find_windows_with_gap
# ---------------------------------------------------------------------------


class TestFindWindowsWithGap:
    def test_no_gaps(self) -> None:
        ref = _nas('ATGATG')
        seq = _nas('ATGATG')
        windows = find_windows_with_gap(ref, seq, 30)
        assert windows == []

    def test_single_gap_window(self) -> None:
        ref = _nas('ATG---ATG')
        seq = _nas('ATGCCCATG')
        windows = find_windows_with_gap(ref, seq, 30)
        assert len(windows) == 1
        assert windows[0] == slice(3, 6)

    def test_two_windows_beyond_distance(self) -> None:
        ref = _nas('ATGATGATGATGATG')
        seq = _nas('-TG--GATGATGATG')
        windows = find_windows_with_gap(ref, seq, 1)
        assert len(windows) == 2

    def test_gaps_merged_within_distance(self) -> None:
        ref = _nas('ATG---ATG---ATG')
        seq = _nas('ATGCCCATGCCCATG')
        windows = find_windows_with_gap(ref, seq, 30)
        assert len(windows) == 1


# ---------------------------------------------------------------------------
# center_expand_positions
# ---------------------------------------------------------------------------


class TestCenterExpandPositions:
    def test_basic_expansion(self) -> None:
        positions = center_expand_positions(6, 0, 12, 3)
        assert positions[0] == 6
        assert len(positions) == len(set(positions))

    def test_scanstart_limit(self) -> None:
        positions = center_expand_positions(0, 3, 12, 3)
        assert all(p >= 3 for p in positions)

    def test_negative_center(self) -> None:
        positions = center_expand_positions(-1, 0, 12, 3)
        assert len(positions) > 0
        assert positions[0] == 0

    def test_center_beyond_range(self) -> None:
        positions = center_expand_positions(100, 0, 12, 3)
        # Should still produce positions covering the valid range
        assert len(positions) > 0
        # All positions should be >= scanstart (0)
        assert all(p >= 0 for p in positions)

    def test_zero_length(self) -> None:
        positions = center_expand_positions(0, 0, 0, 3)
        assert 0 in positions


# ---------------------------------------------------------------------------
# extend_codons_until_gap
# ---------------------------------------------------------------------------


class TestExtendCodonsUntilGap:
    def test_right_no_gaps(self) -> None:
        c1 = [_na('A', 1), _na('T', 2), _na('G', 3)]
        c2 = [_na('C', 4), _na('C', 5), _na('A', 6)]
        _r, _s, count = extend_codons_until_gap([c1, c2], [c1, c2], RIGHT)
        assert count == 2

    def test_right_gap_at_start(self) -> None:
        c1 = [_na('-'), _na('T', 1), _na('G', 2)]
        c2 = [_na('A', 3), _na('T', 4), _na('G', 5)]
        _r, _s, count = extend_codons_until_gap([c1, c2], [c1, c2], RIGHT)
        assert count == 0

    def test_left_direction(self) -> None:
        c1 = [_na('A', 1), _na('T', 2), _na('G', 3)]
        c2 = [_na('C', 4), _na('C', 5), _na('A', 6)]
        _r, _s, count = extend_codons_until_gap([c1, c2], [c1, c2], LEFT)
        assert count == 2

    def test_empty_input(self) -> None:
        _r, _s, count = extend_codons_until_gap([], [], RIGHT)
        assert count == 0


# ---------------------------------------------------------------------------
# parse_gap_placement_score
# ---------------------------------------------------------------------------


class TestParseGapPlacementScore:
    def test_empty_string(self) -> None:
        result = parse_gap_placement_score('')
        assert result == {REFGAP: {}, SEQGAP: {}}

    def test_single_ins(self) -> None:
        result = parse_gap_placement_score('204ins:-5')
        assert result[REFGAP] == {(204, 0): -5}
        assert result[SEQGAP] == {}

    def test_single_del(self) -> None:
        result = parse_gap_placement_score('100del:10')
        assert result[SEQGAP] == {(100, 0): 10}
        assert result[REFGAP] == {}

    def test_with_size(self) -> None:
        result = parse_gap_placement_score('2041/12del:10')
        assert result[SEQGAP] == {(2041, 12): 10}

    def test_multiple_scores(self) -> None:
        result = parse_gap_placement_score('204ins:-5,2041/12del:10')
        assert result[REFGAP] == {(204, 0): -5}
        assert result[SEQGAP] == {(2041, 12): 10}

    def test_invalid_format(self) -> None:
        with pytest.raises(ValueError):
            parse_gap_placement_score('not_valid')

    def test_negative_score(self) -> None:
        result = parse_gap_placement_score('100del:-50')
        assert result[SEQGAP] == {(100, 0): -50}

    def test_zero_score(self) -> None:
        result = parse_gap_placement_score('100del:0')
        assert result[SEQGAP] == {(100, 0): 0}

    def test_trailing_comma(self) -> None:
        result = parse_gap_placement_score('100del:5,')
        assert result[SEQGAP] == {(100, 0): 5}


# ---------------------------------------------------------------------------
# calc_match_score_precomputed
# ---------------------------------------------------------------------------


class TestCalcMatchScorePrecomputed:
    def test_identical_sequences_positive_score(self) -> None:
        nas = _nas('ATGATG')
        other_aas = translate_codons(nas)
        score = calc_match_score_precomputed(nas, nas, other_aas, 0.0)
        assert score > 0

    def test_base_score_applied(self) -> None:
        nas = _nas('ATGATG')
        other_aas = translate_codons(nas)
        score_zero = calc_match_score_precomputed(nas, nas, other_aas, 0.0)
        score_neg = calc_match_score_precomputed(nas, nas, other_aas, -10.0)
        assert score_neg == pytest.approx(score_zero - 10.0)

    def test_mismatched_lower_score(self) -> None:
        ref = _nas('ATGATG')
        seq = _nas('CCCGGG')
        ref_aas = translate_codons(ref)
        score_match = calc_match_score_precomputed(ref, ref, ref_aas, 0.0)
        score_mismatch = calc_match_score_precomputed(seq, ref, ref_aas, 0.0)
        assert score_mismatch < score_match


# ---------------------------------------------------------------------------
# gather_gaps
# ---------------------------------------------------------------------------


class TestGatherGaps:
    def test_no_gaps(self) -> None:
        ref = _nas('ATGATG')
        seq = _nas('ATGATG')
        r, s = gather_gaps(ref, seq, 30)
        assert NAPosition.as_str(r) == 'ATGATG'
        assert NAPosition.as_str(s) == 'ATGATG'

    def test_single_gap(self) -> None:
        ref = _nas('ATG---ATG')
        seq = _nas('ATGCCCATG')
        r, _s = gather_gaps(ref, seq, 30)
        assert NAPosition.count_gaps(r) == 3

    def test_redundant_gaps_removed(self) -> None:
        ref = _nas('ATG---ATG')
        seq = _nas('ATG---ATG')
        r, s = gather_gaps(ref, seq, 30)
        assert NAPosition.count_gaps(r) == 0
        assert NAPosition.count_gaps(s) == 0


# ---------------------------------------------------------------------------
# realign_gaps
# ---------------------------------------------------------------------------


class TestRealignGaps:
    def test_no_gaps_passthrough(self) -> None:
        ref = _nas('ATGATGATG')
        seq = _nas('ATGATGATG')
        gps: dict[int, dict[tuple[int, int], int]] = {REFGAP: {}, SEQGAP: {}}
        r, s = realign_gaps(ref, seq, 30, 10, gps, False, False)
        assert NAPosition.as_str(r) == 'ATGATGATG'
        assert NAPosition.as_str(s) == 'ATGATGATG'

    def test_with_gap(self) -> None:
        ref = _nas('ATGATGATGATG')
        seq = _nas('ATGATG---ATG')
        gps: dict[int, dict[tuple[int, int], int]] = {REFGAP: {}, SEQGAP: {}}
        r, s = realign_gaps(ref, seq, 30, 10, gps, False, False)
        assert len(r) == len(s)
        r_str = NAPosition.as_str(r)
        s_str = NAPosition.as_str(s)
        assert r_str.replace('-', '') == 'ATGATGATGATG'
        assert s_str.replace('-', '') == 'ATGATGATG'

    def test_is_seq_start_flag(self) -> None:
        ref = _nas('ATGATGATG')
        seq = _nas('---ATGATG')
        gps: dict[int, dict[tuple[int, int], int]] = {REFGAP: {}, SEQGAP: {}}
        _r, s = realign_gaps(ref, seq, 30, 10, gps, True, False)
        s_str = NAPosition.as_str(s)
        assert s_str.startswith('---')

    def test_is_seq_end_flag(self) -> None:
        ref = _nas('ATGATGATG')
        seq = _nas('ATGATG---')
        gps: dict[int, dict[tuple[int, int], int]] = {REFGAP: {}, SEQGAP: {}}
        _r, s = realign_gaps(ref, seq, 30, 10, gps, False, True)
        s_str = NAPosition.as_str(s)
        assert s_str.endswith('---')


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------


class TestConstants:
    def test_gap_type_values(self) -> None:
        assert NOGAP == 0
        assert REFGAP == 1
        assert SEQGAP == 2

    def test_direction_values(self) -> None:
        assert LEFT == 0
        assert RIGHT == 1


# ---------------------------------------------------------------------------
# gap_placement_score_callback
# ---------------------------------------------------------------------------


class TestGapPlacementScoreCallback:
    def _make_param(self) -> click.Option:
        return click.Option(['--gap-placement-score'])

    def test_valid_input(self) -> None:
        ctx = click.Context(click.Command('test'))
        param = self._make_param()
        result = gap_placement_score_callback(ctx, param, ('204ins:-5',))
        assert result[REFGAP] == {(204, 0): -5}

    def test_multiple_values_joined(self) -> None:
        ctx = click.Context(click.Command('test'))
        param = self._make_param()
        result = gap_placement_score_callback(ctx, param, ('204ins:-5', '100del:10'))
        assert result[REFGAP] == {(204, 0): -5}
        assert result[SEQGAP] == {(100, 0): 10}

    def test_invalid_raises_bad_option(self) -> None:
        ctx = click.Context(click.Command('test'))
        param = self._make_param()
        with pytest.raises(click.BadOptionUsage):
            gap_placement_score_callback(ctx, param, ('INVALID',))

    def test_empty_param_name_raises(self) -> None:
        ctx = click.Context(click.Command('test'))
        param = self._make_param()
        # Override name to empty string (falsy)
        param.name = ''
        with pytest.raises(click.BadParameter):
            gap_placement_score_callback(ctx, param, ('204ins:-5',))


# ---------------------------------------------------------------------------
# CLI codon_alignment command
# ---------------------------------------------------------------------------


class TestCLICommand:
    def test_command_exists(self) -> None:
        """codon-alignment is registered on the root CLI."""
        from postalign.cli import cli as root_cli

        assert 'codon-alignment' in root_cli.commands

    def test_codon_alignment_returns_processor(self) -> None:
        """Invoke the Click command directly to get a Processor."""
        from postalign.cli import cli as root_cli

        cmd = root_cli.commands['codon-alignment']
        ctx = click.Context(cmd)
        # Invoke with default params using standalone_mode=False
        proc = ctx.invoke(
            cmd,
            min_gap_distance=30,
            window_size=10,
            gap_placement_score={REFGAP: {}, SEQGAP: {}},
            backend='python',
            ref_start=1,
            ref_end=100,
        )
        assert callable(proc)

    def test_ref_start_zero_raises(self) -> None:
        """ref_start < 1 should raise ClickException."""
        from postalign.cli import cli as root_cli

        cmd = root_cli.commands['codon-alignment']
        ctx = click.Context(cmd)
        with pytest.raises(click.ClickException):
            ctx.invoke(
                cmd,
                min_gap_distance=30,
                window_size=10,
                gap_placement_score={REFGAP: {}, SEQGAP: {}},
                backend='python',
                ref_start=0,
                ref_end=100,
            )

    def test_ref_end_too_close_raises(self) -> None:
        """ref_end - 2 < ref_start should raise."""
        from postalign.cli import cli as root_cli

        cmd = root_cli.commands['codon-alignment']
        ctx = click.Context(cmd)
        with pytest.raises(click.ClickException):
            ctx.invoke(
                cmd,
                min_gap_distance=30,
                window_size=10,
                gap_placement_score={REFGAP: {}, SEQGAP: {}},
                backend='python',
                ref_start=10,
                ref_end=10,
            )
