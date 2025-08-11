"""Tests for codon alignment utilities."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
import typer

from postalign.models import NAPosition, Sequence
from postalign.processors.codon_alignment import (
    LEFT,
    REFGAP,
    SEQGAP,
    codon_alignment,
    codon_pairs_group_key,
    extend_codons_until_gap,
    find_best_matches,
    find_first_gap,
    move_gaps_to_center,
    parse_gap_placement_score,
    gap_placement_score_callback,
    calc_match_score,
    separate_gaps_from_nas,
    paired_find_best_matches,
    gather_gaps,
    adjust_gap_placement,
    realign_gaps,
    codon_align,
)


def test_parse_gap_placement_score_valid() -> None:
    """Valid gap placement strings should map to score dictionaries."""

    scores = parse_gap_placement_score("204ins:-5,2041/12del:10")
    assert scores[REFGAP][(204, 0)] == -5
    assert scores[SEQGAP][(2041, 12)] == 10


def test_parse_gap_placement_score_invalid() -> None:
    """Invalid score strings should raise :class:`ValueError`."""

    with pytest.raises(ValueError):
        parse_gap_placement_score("204foo")


def test_parse_gap_placement_score_skips_empty() -> None:
    """Empty entries should be ignored when parsing placement scores."""

    scores = parse_gap_placement_score("204ins:-5,,205del:3,")
    assert scores[REFGAP][(204, 0)] == -5
    assert scores[SEQGAP][(205, 0)] == 3


def test_gap_placement_score_callback_missing_name() -> None:
    """Missing parameter name should raise :class:`BadParameter`."""

    ctx = MagicMock()
    param = MagicMock()
    param.name = None
    with pytest.raises(typer.BadParameter):
        gap_placement_score_callback(ctx, param, ())


def test_codon_alignment_invalid_ref_start() -> None:
    """Invalid reference start should raise :class:`BadParameter`."""

    with pytest.raises(typer.BadParameter):
        codon_alignment(ref_start=0, ref_end=-1)


def test_codon_alignment_invalid_ref_end() -> None:
    """Too-small ref_end should raise :class:`BadParameter`."""

    with pytest.raises(typer.BadParameter):
        codon_alignment(ref_start=5, ref_end=6)


def test_extend_codons_until_gap_right() -> None:
    """Scanning forward should stop before codons containing gaps."""
    from postalign.models import NAPosition
    from postalign.processors.codon_alignment import (
        extend_codons_until_gap,
        RIGHT,
    )

    codon1 = NAPosition.init_from_bytes(b"ATG")
    codon2 = NAPosition.init_from_bytes(b"A-G")
    ref_cd, seq_cd, length = extend_codons_until_gap(
        [codon1, codon2], [codon1, codon2], RIGHT
    )
    assert length == 1
    assert [NAPosition.as_str(c) for c in ref_cd] == ["ATG"]


def test_find_windows_with_gap() -> None:
    """Gap windows should be separated when distance exceeds threshold."""
    from postalign.models import NAPosition
    from postalign.processors.codon_alignment import find_windows_with_gap

    ref = NAPosition.init_from_bytes(b"AAA---BBB---CCC")
    seq = NAPosition.init_from_bytes(b"AAA---BBB---CCC")
    windows = find_windows_with_gap(ref, seq, 1)
    assert windows == [slice(3, 6), slice(9, 12)]


def test_find_first_gap_none() -> None:
    """Sequences without gaps should return ``-1`` for first gap."""
    from postalign.models import NAPosition
    from postalign.processors.codon_alignment import find_first_gap

    nas = NAPosition.init_from_bytes(b"ATG")
    assert find_first_gap(nas) == -1


def test_move_gap_to_codon_end() -> None:
    """Gaps within codons should be shifted to the end."""
    from postalign.models import NAPosition
    from postalign.processors.codon_alignment import move_gap_to_codon_end

    codons = [NAPosition.init_from_bytes(b"A-A")]
    moved = move_gap_to_codon_end(codons)
    assert [NAPosition.as_str(c) for c in moved] == ["AA-"]


def test_remove_redundant_gaps() -> None:
    """Matched gaps across sequences should be pruned."""
    from postalign.models import NAPosition
    from postalign.processors.codon_alignment import remove_redundant_gaps

    ref = NAPosition.init_from_bytes(b"A--B")
    seq = NAPosition.init_from_bytes(b"A--B")
    new_ref, new_seq = remove_redundant_gaps(ref, seq)
    assert NAPosition.as_str(new_ref) == "AB"
    assert NAPosition.as_str(new_seq) == "AB"


def test_extend_codons_until_gap_left() -> None:
    """Scanning backward should stop before codons containing gaps."""

    codon1 = NAPosition.init_from_bytes(b"A-G")
    codon2 = NAPosition.init_from_bytes(b"ATG")
    ref_cd, seq_cd, length = extend_codons_until_gap(
        [codon1, codon2], [codon1, codon2], LEFT
    )
    assert length == 1
    assert [NAPosition.as_str(c) for c in ref_cd] == ["ATG"]


def test_find_first_gap_found() -> None:
    """Sequences with gaps should return the first gap index."""

    nas = NAPosition.init_from_bytes(b"AT-G")
    assert find_first_gap(nas) == 2


def test_calc_match_score() -> None:
    """Match scoring should include nucleotide and amino-acid contributions."""

    mynas = NAPosition.init_from_bytes(b"ATG")
    othernas = NAPosition.init_from_bytes(b"ATG")
    assert calc_match_score(mynas, othernas, 0.0) == pytest.approx(8.0)


def test_separate_gaps_from_nas() -> None:
    """Separating gaps should partition sequences into nongaps and gaps."""

    nas = NAPosition.init_from_bytes(b"A-G")
    nongaps, gaps = separate_gaps_from_nas(nas)
    assert NAPosition.as_str(nongaps) == "AG"
    assert NAPosition.as_str(gaps) == "-"


def test_find_best_matches_fallback() -> None:
    """If no placement improves score, original nucleotides are returned."""

    mynas = NAPosition.init_from_bytes(b"A-")
    othernas = NAPosition.init_from_bytes(b"A")
    result = find_best_matches(
        mynas, othernas, set(), REFGAP, {}, False, False
    )
    assert NAPosition.as_str(result) == "A"


def test_find_best_matches_prefers_scored_position() -> None:
    """Gap placement scores should influence the chosen insertion index."""

    mynas = NAPosition.init_from_bytes(b"-AAAA")
    othernas = NAPosition.init_from_bytes(b"AAAA")
    gps = {(4, 1): 5}
    with patch(
        "postalign.processors.codon_alignment.calc_match_score",
        return_value=0.0,
    ):
        result = find_best_matches(
            mynas, othernas, set(), SEQGAP, gps, False, False
        )
    assert NAPosition.as_str(result) == "AAA-A"


def test_find_best_matches_any_size_score_and_start() -> None:
    """Generic position scores should apply when size-specific score absent."""

    mynas = NAPosition.init_from_bytes(b"AAAA-")
    othernas = NAPosition.init_from_bytes(b"AAAA")
    gps = {(1, 0): 3}
    with patch(
        "postalign.processors.codon_alignment.calc_match_score",
        return_value=0.0,
    ):
        result = find_best_matches(
            mynas, othernas, set(), SEQGAP, gps, True, False
        )
    assert NAPosition.as_str(result) == "-AAAA"


def test_paired_find_best_matches_refgap() -> None:
    """REFGAP type should delegate to ``find_best_matches`` with ref NAs."""

    ref = NAPosition.init_from_bytes(b"A-A")
    seq = NAPosition.init_from_bytes(b"AAA")
    gps: dict[int, dict[tuple[int, int], int]] = {REFGAP: {}, SEQGAP: {}}
    with patch(
        "postalign.processors.codon_alignment.find_best_matches"
    ) as fbm:
        fbm.return_value = ref
        paired_find_best_matches(ref, seq, REFGAP, gps, False, False)
        fbm.assert_called_once()


def test_paired_find_best_matches_seqgap() -> None:
    """SEQGAP type should delegate using sequence nucleotides."""

    ref = NAPosition.init_from_bytes(b"AAA")
    seq = NAPosition.init_from_bytes(b"A-A")
    gps: dict[int, dict[tuple[int, int], int]] = {REFGAP: {}, SEQGAP: {}}
    with patch(
        "postalign.processors.codon_alignment.find_best_matches"
    ) as fbm:
        fbm.return_value = seq
        paired_find_best_matches(ref, seq, SEQGAP, gps, True, True)
        fbm.assert_called_once()


def test_codon_pairs_group_key() -> None:
    """Codon pairs should classify gap types correctly."""

    refgap_cd = NAPosition.init_from_bytes(b"A--")
    seq_cd = NAPosition.init_from_bytes(b"AAA")
    assert codon_pairs_group_key((0, (refgap_cd, seq_cd))) == REFGAP

    ref_cd2 = NAPosition.init_from_bytes(b"AAA")
    seqgap_cd = NAPosition.init_from_bytes(b"A--")
    assert codon_pairs_group_key((0, (ref_cd2, seqgap_cd))) == SEQGAP

    ref_cd3 = NAPosition.init_from_bytes(b"AAA")
    seq_cd3 = NAPosition.init_from_bytes(b"AAA")
    assert codon_pairs_group_key((0, (ref_cd3, seq_cd3))) == 0


def test_move_gaps_to_center() -> None:
    """Gaps should be relocated to the center of the sequence."""

    nas = NAPosition.init_from_bytes(b"--AB")
    centered = move_gaps_to_center(nas)
    assert NAPosition.as_str(centered) == "A--B"


def test_gather_gaps_relocates_gaps() -> None:
    """Gap windows should be centered and redundant gaps removed."""

    ref = NAPosition.init_from_bytes(b"ABC")
    seq = NAPosition.init_from_bytes(b"A--BC")
    out_ref, out_seq = gather_gaps(ref, seq, 0)
    assert NAPosition.as_str(out_ref) == "ABC"
    assert NAPosition.as_str(out_seq) == "A--BC"


def test_adjust_gap_placement_extends_windows() -> None:
    """Codon windows should extend to include flanking codons."""

    from itertools import chain

    refcodons = [
        NAPosition.init_from_bytes(b"AAA"),
        NAPosition.init_from_bytes(b"A--"),
        NAPosition.init_from_bytes(b"CCC"),
    ]
    seqcodons = [
        NAPosition.init_from_bytes(b"AAA"),
        NAPosition.init_from_bytes(b"AAA"),
        NAPosition.init_from_bytes(b"CCC"),
    ]
    gps: dict[int, dict[tuple[int, int], int]] = {REFGAP: {}, SEQGAP: {}}
    with patch(
        "postalign.processors.codon_alignment.extend_codons_until_gap",
        side_effect=[
            (
                [NAPosition.init_from_bytes(b"GGG")],
                [NAPosition.init_from_bytes(b"GGG")],
                1,
            ),
            (
                [NAPosition.init_from_bytes(b"TTT")],
                [NAPosition.init_from_bytes(b"TTT")],
                1,
            ),
        ],
    ) as ext, patch(
        "postalign.processors.codon_alignment.paired_find_best_matches",
        return_value=(
            list(chain(*refcodons)),
            list(chain(*seqcodons)),
        ),
    ) as pfbm:
        adjust_gap_placement(refcodons, seqcodons, 1, gps, False, False)
    assert ext.call_count == 2
    pfbm.assert_called_once()


def test_realign_gaps_invokes_helpers() -> None:
    """Realignment should delegate to gather and adjustment helpers."""

    ref = NAPosition.init_from_bytes(b"A--A")
    seq = NAPosition.init_from_bytes(b"AA--")
    gps: dict[int, dict[tuple[int, int], int]] = {REFGAP: {}, SEQGAP: {}}
    with patch(
        "postalign.processors.codon_alignment.gather_gaps",
        return_value=(ref, seq),
    ) as gg, patch(
        "postalign.processors.codon_alignment.group_by_codons",
        return_value=([ref], [seq]),
    ) as gbc, patch(
        "postalign.processors.codon_alignment.adjust_gap_placement",
        return_value=([ref], [seq]),
    ) as agp, patch(
        "postalign.processors.codon_alignment.move_gap_to_codon_end",
        return_value=[ref],
    ) as mg:
        realign_gaps(ref, seq, 1, 1, gps, True, True)
    gg.assert_called_once()
    gbc.assert_called_once()
    agp.assert_called_once()
    mg.assert_called_once()


def test_codon_align_no_gaps_returns_input() -> None:
    """Inputs without gaps should pass through codon_align unchanged."""

    ref = Sequence(
        header="r",
        description="",
        seqtext=NAPosition.init_from_bytes(b"AAA"),
        seqid=1,
        seqtype=NAPosition,
        abs_seqstart=0,
    )
    seq = Sequence(
        header="s",
        description="",
        seqtext=NAPosition.init_from_bytes(b"AAA"),
        seqid=1,
        seqtype=NAPosition,
        abs_seqstart=0,
    )
    out_ref, out_seq = codon_align(
        ref,
        seq,
        1,
        1,
        {REFGAP: {}, SEQGAP: {}},
        1,
        3,
    )
    assert out_ref.seqtext_as_str == "AAA"
    assert out_seq.seqtext_as_str == "AAA"


def test_codon_alignment_processor_dispatches() -> None:
    """Processor should skip empty sequences and call ``codon_align``."""

    proc = codon_alignment()
    empty_ref = Sequence(
        header="r1",
        description="",
        seqtext=[],
        seqid=1,
        seqtype=NAPosition,
        abs_seqstart=0,
    )
    empty_seq = Sequence(
        header="s1",
        description="",
        seqtext=[],
        seqid=1,
        seqtype=NAPosition,
        abs_seqstart=0,
    )
    full_ref = Sequence(
        header="r2",
        description="",
        seqtext=NAPosition.init_from_bytes(b"AAA"),
        seqid=2,
        seqtype=NAPosition,
        abs_seqstart=0,
    )
    full_seq = Sequence(
        header="s2",
        description="",
        seqtext=NAPosition.init_from_bytes(b"AAA"),
        seqid=2,
        seqtype=NAPosition,
        abs_seqstart=0,
    )
    with patch(
        "postalign.processors.codon_alignment.codon_align",
        return_value=(full_ref, full_seq),
    ) as ca:
        result = list(
            proc(
                iter([(empty_ref, empty_seq), (full_ref, full_seq)]),
                [],
            )
        )
    ca.assert_called_once()
    assert result[0] == (empty_ref, empty_seq)
    assert result[1] == (full_ref, full_seq)


def test_gap_placement_score_callback_invalid_value() -> None:
    """Invalid callback values should raise :class:`BadParameter`."""

    ctx = MagicMock()
    param = MagicMock()
    param.name = "gps"
    with pytest.raises(typer.BadParameter):
        gap_placement_score_callback(ctx, param, ("204foo",))


def test_gap_placement_score_callback_str_value() -> None:
    """String values should be accepted when Typer lacks ``multiple``."""

    ctx = MagicMock()
    param = MagicMock()
    param.name = "gap_placement_score"
    result = gap_placement_score_callback(
        ctx, param, "204ins:-5,2041/12del:10"
    )
    assert result[REFGAP][(204, 0)] == -5
    assert result[SEQGAP][(2041, 12)] == 10
