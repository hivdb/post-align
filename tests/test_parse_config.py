"""Tests for parse_gap_placement_score() — stable public API."""

import pytest

from postalign.processors.codon_alignment import (
    REFGAP,
    SEQGAP,
    parse_gap_placement_score,
)


class TestParseGapPlacementScore:

    def test_empty_string(self) -> None:
        result = parse_gap_placement_score('')
        assert result == {REFGAP: {}, SEQGAP: {}}

    def test_single_insertion_score(self) -> None:
        result = parse_gap_placement_score('204ins:-5')
        assert result[REFGAP] == {(204, 0): -5}
        assert result[SEQGAP] == {}

    def test_single_deletion_score(self) -> None:
        result = parse_gap_placement_score('2041/12del:10')
        assert result[SEQGAP] == {(2041, 12): 10}
        assert result[REFGAP] == {}

    def test_multiple_scores(self) -> None:
        result = parse_gap_placement_score('204ins:-5,2041/12del:10')
        assert result[REFGAP] == {(204, 0): -5}
        assert result[SEQGAP] == {(2041, 12): 10}

    def test_positive_insertion_score(self) -> None:
        result = parse_gap_placement_score('100ins:20')
        assert result[REFGAP] == {(100, 0): 20}

    def test_deletion_with_size(self) -> None:
        result = parse_gap_placement_score('500/6del:-3')
        assert result[SEQGAP] == {(500, 6): -3}

    def test_insertion_with_size(self) -> None:
        result = parse_gap_placement_score('300/9ins:15')
        assert result[REFGAP] == {(300, 9): 15}

    def test_multiple_same_type(self) -> None:
        result = parse_gap_placement_score('100ins:5,200ins:-3')
        assert result[REFGAP] == {(100, 0): 5, (200, 0): -3}

    def test_trailing_comma(self) -> None:
        result = parse_gap_placement_score('204ins:-5,')
        assert result[REFGAP] == {(204, 0): -5}

    def test_invalid_format_raises(self) -> None:
        with pytest.raises(ValueError, match='invalid argument value'):
            parse_gap_placement_score('notvalid')

    def test_invalid_mixed_raises(self) -> None:
        with pytest.raises(ValueError, match='invalid argument value'):
            parse_gap_placement_score('204ins:-5,bad')

    def test_zero_score(self) -> None:
        result = parse_gap_placement_score('100ins:0')
        assert result[REFGAP] == {(100, 0): 0}
