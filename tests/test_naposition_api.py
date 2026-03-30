"""Comprehensive API tests for NAPosition.

Tests every field, constructor parameter, dunder method, and
static/class method on NAPosition to serve as a compatibility
contract when migrating from Cython to Rust #[pyclass].

All tests here must pass *unchanged* after the migration.
"""

from copy import copy

from postalign.models.na_position import (
    NAPosition,
    enumerate_seq_pos,
)
from postalign.models.position_flag import PositionFlag

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _na(
    ch: str,
    pos: int = -1,
    flag: PositionFlag = PositionFlag.NONE,
) -> NAPosition:
    return NAPosition(ord(ch), pos, flag)


def _nas(s: str, start_pos: int = 1) -> list[NAPosition]:
    """Build NAPosition list from string. Gaps get pos=-1."""
    result: list[NAPosition] = []
    pos = start_pos
    for ch in s.upper():
        n = ord(ch)
        if ch in '-.':
            result.append(NAPosition(n, -1, PositionFlag.NONE))
        else:
            result.append(NAPosition(n, pos, PositionFlag.NONE))
            pos += 1
    return result


# ===========================================================================
# 1. Constructor
# ===========================================================================


class TestConstructor:
    def test_basic_construction(self) -> None:
        na = NAPosition(ord('A'), 1, PositionFlag.NONE)
        assert na.notation == ord('A')
        assert na.pos == 1
        assert na.flag == PositionFlag.NONE

    def test_with_payload(self) -> None:
        payload = {'key': 'value'}
        na = NAPosition(ord('A'), 1, PositionFlag.NONE, payload)
        assert na.payload is payload

    def test_default_payload_is_none(self) -> None:
        na = NAPosition(ord('A'), 1, PositionFlag.NONE)
        assert na.payload is None

    def test_gap_construction(self) -> None:
        na = NAPosition(ord('-'), -1, PositionFlag.NONE)
        assert na.is_gap is True
        assert na.pos == -1

    def test_dot_gap_construction(self) -> None:
        na = NAPosition(ord('.'), -1, PositionFlag.NONE)
        assert na.is_gap is True

    def test_nongap_construction(self) -> None:
        na = NAPosition(ord('A'), 1, PositionFlag.NONE)
        assert na.is_gap is False

    def test_flag_as_int(self) -> None:
        na = NAPosition(ord('A'), 1, 0x10)
        assert na.flag == PositionFlag.UNALIGNED
        assert na.flag & PositionFlag.UNALIGNED


# ===========================================================================
# 2. Field access (get/set)
# ===========================================================================


class TestFieldAccess:
    def test_get_notation(self) -> None:
        na = _na('A', 1)
        assert na.notation == ord('A')

    def test_get_pos(self) -> None:
        na = _na('A', 42)
        assert na.pos == 42

    def test_get_flag(self) -> None:
        na = NAPosition(ord('A'), 1, PositionFlag.UNALIGNED)
        assert na.flag == PositionFlag.UNALIGNED

    def test_get_is_gap_true(self) -> None:
        assert _na('-').is_gap is True

    def test_get_is_gap_false(self) -> None:
        assert _na('A', 1).is_gap is False

    def test_get_payload(self) -> None:
        obj = [1, 2, 3]
        na = NAPosition(ord('A'), 1, PositionFlag.NONE, obj)
        assert na.payload is obj


# ===========================================================================
# 3. Dunder methods
# ===========================================================================


class TestDunderMethods:
    def test_str(self) -> None:
        assert str(_na('A', 1)) == 'A'
        assert str(_na('-')) == '-'
        assert str(_na('G', 5)) == 'G'

    def test_bytes(self) -> None:
        assert bytes(_na('A', 1)) == b'A'
        assert bytes(_na('-')) == b'-'

    def test_repr(self) -> None:
        na = NAPosition(ord('A'), 1, PositionFlag.NONE)
        r = repr(na)
        assert 'NAPosition' in r
        assert str(ord('A')) in r

    def test_repr_with_payload(self) -> None:
        na = NAPosition(ord('A'), 1, PositionFlag.NONE, 'hello')
        r = repr(na)
        assert 'hello' in r

    def test_copy(self) -> None:
        original = NAPosition(ord('A'), 1, PositionFlag.UNALIGNED, [1, 2])
        copied = copy(original)
        assert copied.notation == original.notation
        assert copied.pos == original.pos
        assert copied.flag == original.flag
        assert copied.is_gap == original.is_gap
        assert copied.payload == original.payload
        # Should be a different object
        assert copied is not original


# ===========================================================================
# 4. PositionFlag interop
# ===========================================================================


class TestPositionFlagInterop:
    def test_flag_bitwise_or(self) -> None:
        na = _na('A', 1)
        na.flag |= PositionFlag.UNALIGNED
        assert na.flag & PositionFlag.UNALIGNED

    def test_flag_bitwise_or_multiple(self) -> None:
        na = _na('A', 1)
        na.flag |= PositionFlag.UNALIGNED
        na.flag |= PositionFlag.TRIM_BY_SEQ
        assert na.flag & PositionFlag.UNALIGNED
        assert na.flag & PositionFlag.TRIM_BY_SEQ

    def test_flag_starts_as_none(self) -> None:
        na = _na('A', 1)
        assert na.flag == PositionFlag.NONE
        assert not (na.flag & PositionFlag.UNALIGNED)

    def test_flag_as_int_comparison(self) -> None:
        na = NAPosition(ord('A'), 1, 0x11)
        assert na.flag & 0x10  # UNALIGNED
        assert na.flag & 0x01  # TRIM_BY_SEQ


# ===========================================================================
# 5. isinstance checks
# ===========================================================================


class TestIsInstance:
    def test_isinstance_naposition(self) -> None:
        na = _na('A', 1)
        assert isinstance(na, NAPosition)

    def test_type_is_naposition(self) -> None:
        na = _na('A', 1)
        assert type(na) is NAPosition


# ===========================================================================
# 6. Class method: init_gaps
# ===========================================================================


class TestInitGaps:
    def test_basic(self) -> None:
        gaps = NAPosition.init_gaps(3)
        assert len(gaps) == 3
        for g in gaps:
            assert g.is_gap is True
            assert g.pos == -1
            assert g.notation == ord('-')
            assert g.flag == PositionFlag.NONE

    def test_zero_length(self) -> None:
        gaps = NAPosition.init_gaps(0)
        assert gaps == []

    def test_single(self) -> None:
        gaps = NAPosition.init_gaps(1)
        assert len(gaps) == 1
        assert gaps[0].is_gap is True


# ===========================================================================
# 7. Class method: init_from_bytes
# ===========================================================================


class TestInitFromBytes:
    def test_basic(self) -> None:
        nas = NAPosition.init_from_bytes(b'ATG')
        assert len(nas) == 3
        assert nas[0].notation == ord('A')
        assert nas[0].pos == 1
        assert nas[1].notation == ord('T')
        assert nas[1].pos == 2
        assert nas[2].notation == ord('G')
        assert nas[2].pos == 3

    def test_uppercase_conversion(self) -> None:
        nas = NAPosition.init_from_bytes(b'atg')
        assert nas[0].notation == ord('A')
        assert nas[1].notation == ord('T')
        assert nas[2].notation == ord('G')

    def test_with_gaps(self) -> None:
        nas = NAPosition.init_from_bytes(b'A-G')
        assert len(nas) == 3
        assert nas[0].pos == 1
        assert nas[1].pos == -1
        assert nas[1].is_gap is True
        assert nas[2].pos == 2

    def test_with_payload(self) -> None:
        nas = NAPosition.init_from_bytes(b'ATG', ['p1', 'p2', 'p3'])
        assert nas[0].payload == 'p1'
        assert nas[1].payload == 'p2'
        assert nas[2].payload == 'p3'

    def test_payload_none_default(self) -> None:
        nas = NAPosition.init_from_bytes(b'ATG')
        for na in nas:
            assert na.payload is None

    def test_partial_payload(self) -> None:
        nas = NAPosition.init_from_bytes(b'ATG', ['p1'])
        assert nas[0].payload == 'p1'
        assert nas[1].payload is None
        assert nas[2].payload is None

    def test_empty(self) -> None:
        nas = NAPosition.init_from_bytes(b'')
        assert nas == []

    def test_bytearray_input(self) -> None:
        nas = NAPosition.init_from_bytes(bytearray(b'ATG'))
        assert len(nas) == 3
        assert nas[0].notation == ord('A')


# ===========================================================================
# 8. Static method: min_pos / max_pos
# ===========================================================================


class TestMinMaxPos:
    def test_min_pos_basic(self) -> None:
        nas = _nas('ATG')
        assert NAPosition.min_pos(nas) == 1

    def test_max_pos_basic(self) -> None:
        nas = _nas('ATG')
        assert NAPosition.max_pos(nas) == 3

    def test_min_pos_with_leading_gaps(self) -> None:
        nas = _nas('--ATG')
        assert NAPosition.min_pos(nas) == 1

    def test_max_pos_with_trailing_gaps(self) -> None:
        nas = _nas('ATG--')
        assert NAPosition.max_pos(nas) == 3

    def test_min_pos_all_gaps(self) -> None:
        nas = _nas('---')
        assert NAPosition.min_pos(nas) == -1

    def test_max_pos_all_gaps(self) -> None:
        nas = _nas('---')
        assert NAPosition.max_pos(nas) == -1

    def test_min_pos_empty(self) -> None:
        assert NAPosition.min_pos([]) == -1

    def test_max_pos_empty(self) -> None:
        assert NAPosition.max_pos([]) == -1

    def test_min_pos_single_nongap(self) -> None:
        nas = _nas('--A--')
        assert NAPosition.min_pos(nas) == 1

    def test_max_pos_single_nongap(self) -> None:
        nas = _nas('--A--')
        assert NAPosition.max_pos(nas) == 1


# ===========================================================================
# 9. Static method: min_nongap_index / max_nongap_index
# ===========================================================================


class TestNongapIndex:
    def test_min_nongap_index_basic(self) -> None:
        nas = _nas('ATG')
        assert NAPosition.min_nongap_index(nas) == 0

    def test_min_nongap_index_with_leading_gaps(self) -> None:
        nas = _nas('--ATG')
        assert NAPosition.min_nongap_index(nas) == 2

    def test_min_nongap_index_all_gaps(self) -> None:
        nas = _nas('---')
        assert NAPosition.min_nongap_index(nas) == -1

    def test_min_nongap_index_with_start(self) -> None:
        nas = _nas('A-TG')
        assert NAPosition.min_nongap_index(nas, start=2) == 2

    def test_min_nongap_index_with_start_and_stop(self) -> None:
        nas = _nas('A-TG')
        assert NAPosition.min_nongap_index(nas, start=1, stop=2) == -1

    def test_max_nongap_index_basic(self) -> None:
        nas = _nas('ATG')
        assert NAPosition.max_nongap_index(nas) == 2

    def test_max_nongap_index_with_trailing_gaps(self) -> None:
        nas = _nas('ATG--')
        assert NAPosition.max_nongap_index(nas) == 2

    def test_max_nongap_index_all_gaps(self) -> None:
        nas = _nas('---')
        assert NAPosition.max_nongap_index(nas) == -1

    def test_max_nongap_index_with_start_stop(self) -> None:
        nas = _nas('ATG--ATG')
        assert NAPosition.max_nongap_index(nas, start=0, stop=3) == 2

    def test_min_nongap_index_empty(self) -> None:
        assert NAPosition.min_nongap_index([]) == -1

    def test_max_nongap_index_empty(self) -> None:
        assert NAPosition.max_nongap_index([]) == -1


# ===========================================================================
# 10. Static method: count_gaps / count_nongaps
# ===========================================================================


class TestCountGaps:
    def test_count_gaps_no_gaps(self) -> None:
        assert NAPosition.count_gaps(_nas('ATG')) == 0

    def test_count_gaps_all_gaps(self) -> None:
        assert NAPosition.count_gaps(_nas('---')) == 3

    def test_count_gaps_mixed(self) -> None:
        assert NAPosition.count_gaps(_nas('A-G')) == 1

    def test_count_gaps_empty(self) -> None:
        assert NAPosition.count_gaps([]) == 0

    def test_count_nongaps_no_gaps(self) -> None:
        assert NAPosition.count_nongaps(_nas('ATG')) == 3

    def test_count_nongaps_all_gaps(self) -> None:
        assert NAPosition.count_nongaps(_nas('---')) == 0

    def test_count_nongaps_mixed(self) -> None:
        assert NAPosition.count_nongaps(_nas('A-G')) == 2

    def test_count_nongaps_empty(self) -> None:
        assert NAPosition.count_nongaps([]) == 0


# ===========================================================================
# 11. Static method: any_has_gap / all_have_gap
# ===========================================================================


class TestHasGap:
    def test_any_has_gap_true(self) -> None:
        assert NAPosition.any_has_gap(_nas('A-G')) is True

    def test_any_has_gap_false(self) -> None:
        assert NAPosition.any_has_gap(_nas('ATG')) is False

    def test_any_has_gap_empty(self) -> None:
        assert NAPosition.any_has_gap([]) is False

    def test_all_have_gap_true(self) -> None:
        assert NAPosition.all_have_gap(_nas('---')) is True

    def test_all_have_gap_false(self) -> None:
        assert NAPosition.all_have_gap(_nas('A--')) is False

    def test_all_have_gap_empty(self) -> None:
        assert NAPosition.all_have_gap([]) is True

    def test_dot_gap_detected(self) -> None:
        assert NAPosition.any_has_gap(_nas('A.G')) is True


# ===========================================================================
# 12. Static method: set_flag / any_has_flag / all_have_flag
# ===========================================================================


class TestFlagMethods:
    def test_set_flag(self) -> None:
        nas = _nas('ATG')
        NAPosition.set_flag(nas, PositionFlag.UNALIGNED)
        for na in nas:
            assert na.flag & PositionFlag.UNALIGNED

    def test_set_flag_preserves_existing(self) -> None:
        nas = _nas('ATG')
        NAPosition.set_flag(nas, PositionFlag.TRIM_BY_SEQ)
        NAPosition.set_flag(nas, PositionFlag.UNALIGNED)
        for na in nas:
            assert na.flag & PositionFlag.TRIM_BY_SEQ
            assert na.flag & PositionFlag.UNALIGNED

    def test_any_has_flag_true(self) -> None:
        nas = _nas('ATG')
        nas[1].flag |= PositionFlag.UNALIGNED
        assert NAPosition.any_has_flag(nas, PositionFlag.UNALIGNED) is True

    def test_any_has_flag_false(self) -> None:
        nas = _nas('ATG')
        assert NAPosition.any_has_flag(nas, PositionFlag.UNALIGNED) is False

    def test_all_have_flag_true(self) -> None:
        nas = _nas('ATG')
        NAPosition.set_flag(nas, PositionFlag.UNALIGNED)
        assert NAPosition.all_have_flag(nas, PositionFlag.UNALIGNED) is True

    def test_all_have_flag_false(self) -> None:
        nas = _nas('ATG')
        nas[0].flag |= PositionFlag.UNALIGNED
        assert NAPosition.all_have_flag(nas, PositionFlag.UNALIGNED) is False

    def test_any_has_flag_empty(self) -> None:
        assert NAPosition.any_has_flag([], PositionFlag.UNALIGNED) is False

    def test_all_have_flag_empty(self) -> None:
        assert NAPosition.all_have_flag([], PositionFlag.UNALIGNED) is True


# ===========================================================================
# 13. Static method: remove_gaps
# ===========================================================================


class TestRemoveGaps:
    def test_remove_gaps_basic(self) -> None:
        nas = _nas('A-G')
        result = NAPosition.remove_gaps(nas)
        assert len(result) == 2
        assert all(not na.is_gap for na in result)

    def test_remove_gaps_no_gaps(self) -> None:
        nas = _nas('ATG')
        result = NAPosition.remove_gaps(nas)
        assert len(result) == 3

    def test_remove_gaps_all_gaps(self) -> None:
        nas = _nas('---')
        result = NAPosition.remove_gaps(nas)
        assert len(result) == 0

    def test_remove_gaps_empty(self) -> None:
        assert NAPosition.remove_gaps([]) == []


# ===========================================================================
# 14. Static method: as_bytes / as_str
# ===========================================================================


class TestAsBytes:
    def test_as_bytes_basic(self) -> None:
        nas = _nas('ATG')
        assert NAPosition.as_bytes(nas) == b'ATG'

    def test_as_bytes_with_gaps(self) -> None:
        nas = _nas('A-G')
        assert NAPosition.as_bytes(nas) == b'A-G'

    def test_as_bytes_empty(self) -> None:
        assert NAPosition.as_bytes([]) == b''

    def test_as_str_basic(self) -> None:
        nas = _nas('ATG')
        assert NAPosition.as_str(nas) == 'ATG'

    def test_as_str_with_gaps(self) -> None:
        nas = _nas('ATG---CCA')
        assert NAPosition.as_str(nas) == 'ATG---CCA'

    def test_as_str_empty(self) -> None:
        assert NAPosition.as_str([]) == ''


# ===========================================================================
# 15. Static method: posrange2indexrange
# ===========================================================================


class TestPosrange2Indexrange:
    def test_basic_range(self) -> None:
        nas = _nas('ATGATG')
        start, end = NAPosition.posrange2indexrange(nas, 1, 6)
        assert start == 0
        assert end == 6

    def test_partial_range(self) -> None:
        nas = _nas('ATGATG')
        start, end = NAPosition.posrange2indexrange(nas, 2, 5)
        assert start == 1
        assert end == 5

    def test_with_leading_gaps(self) -> None:
        nas = _nas('--ATG')
        start, end = NAPosition.posrange2indexrange(nas, 1, 3)
        assert start == 2
        assert end == 5

    def test_include_boundary_gaps(self) -> None:
        nas = _nas('--ATG--')
        start, end = NAPosition.posrange2indexrange(nas, 1, 3, include_boundary_gaps=True)
        assert start == 0  # extended to include leading gaps
        assert end == 7  # extended to include trailing gaps

    def test_no_boundary_gaps(self) -> None:
        nas = _nas('--ATG--')
        start, end = NAPosition.posrange2indexrange(nas, 1, 3, include_boundary_gaps=False)
        assert start == 2
        assert end == 5

    def test_all_gaps(self) -> None:
        nas = _nas('---')
        start, end = NAPosition.posrange2indexrange(nas, 1, 3)
        assert (start, end) == (0, 0)

    def test_empty(self) -> None:
        start, end = NAPosition.posrange2indexrange([], 1, 3)
        assert (start, end) == (0, 0)

    def test_pos_beyond_max(self) -> None:
        nas = _nas('ATG')
        start, end = NAPosition.posrange2indexrange(nas, 10, 20)
        # pos_start > max_pos, returns last index
        assert start == end


# ===========================================================================
# 16. Standalone function: enumerate_seq_pos
# ===========================================================================


class TestEnumerateSeqPos:
    def test_basic(self) -> None:
        result = enumerate_seq_pos(b'ATG')
        assert result == [1, 2, 3]

    def test_with_gaps(self) -> None:
        result = enumerate_seq_pos(b'A-G')
        assert result == [1, -1, 2]

    def test_all_gaps(self) -> None:
        result = enumerate_seq_pos(b'---')
        assert result == [-1, -1, -1]

    def test_dot_gaps(self) -> None:
        result = enumerate_seq_pos(b'A.G')
        assert result == [1, -1, 2]

    def test_empty(self) -> None:
        result = enumerate_seq_pos(b'')
        assert result == []


# ===========================================================================
# 17. Type token usage (seqtype pattern)
# ===========================================================================


class TestTypeToken:
    def test_callable_init_gaps(self) -> None:
        """NAPosition can be used as seqtype.init_gaps()."""
        seqtype = NAPosition
        gaps = seqtype.init_gaps(3)
        assert len(gaps) == 3

    def test_callable_init_from_bytes(self) -> None:
        """NAPosition can be used as seqtype.init_from_bytes()."""
        seqtype = NAPosition
        nas = seqtype.init_from_bytes(b'ATG')
        assert len(nas) == 3

    def test_callable_set_flag(self) -> None:
        """NAPosition can be used as seqtype.set_flag()."""
        seqtype = NAPosition
        nas = _nas('ATG')
        seqtype.set_flag(nas, PositionFlag.UNALIGNED)
        assert all(na.flag & PositionFlag.UNALIGNED for na in nas)

    def test_callable_count_nongaps(self) -> None:
        """NAPosition can be used as seqtype.count_nongaps()."""
        seqtype = NAPosition
        assert seqtype.count_nongaps(_nas('A-G')) == 2

    def test_callable_as_str(self) -> None:
        """NAPosition can be used as seqtype.as_str()."""
        seqtype = NAPosition
        assert seqtype.as_str(_nas('ATG')) == 'ATG'

    def test_callable_count_gaps(self) -> None:
        """NAPosition can be used as seqtype.count_gaps()."""
        seqtype = NAPosition
        assert seqtype.count_gaps(_nas('A-G')) == 1

    def test_callable_all_have_gap(self) -> None:
        """NAPosition can be used as seqtype.all_have_gap()."""
        seqtype = NAPosition
        assert seqtype.all_have_gap(_nas('---')) is True
        assert seqtype.all_have_gap(_nas('A--')) is False
