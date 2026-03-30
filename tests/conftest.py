"""Shared test fixtures and helpers for codon alignment tests."""

from collections.abc import Callable
from typing import Any

from postalign.models.na_position import NAPosition
from postalign.models.position_flag import PositionFlag
from postalign.models.sequence import Sequence
from postalign.processors.codon_alignment import REFGAP, SEQGAP

DEFAULT_GAP_PLACEMENT_SCORE: dict[int, dict[tuple[int, int], int]] = {
    REFGAP: {},
    SEQGAP: {},
}


def make_na_positions(seq_str: str, start_pos: int = 1) -> list[NAPosition]:
    """Build a list of NAPosition from a string like 'ATG---CCA'.

    Non-gap characters get sequential positions starting from start_pos.
    Gap characters ('-' or '.') get pos=-1.
    """
    positions: list[NAPosition] = []
    pos = start_pos
    for ch in seq_str.upper():
        notation = ord(ch)
        if ch in '-.':
            positions.append(NAPosition(notation, -1, PositionFlag.NONE))
        else:
            positions.append(NAPosition(notation, pos, PositionFlag.NONE))
            pos += 1
    return positions


def make_sequence(
    seq_str: str,
    header: str = 'test',
    seqid: int = 0,
    start_pos: int = 1,
) -> Sequence:
    """Build a Sequence[NAPosition] from a plain string like 'ATG---CCA'."""
    nas = make_na_positions(seq_str, start_pos=start_pos)
    return Sequence(
        header=header,
        description='',
        seqtext=nas,
        seqid=seqid,
        seqtype=NAPosition,
        abs_seqstart=0,
        skip_invalid=True,
    )


def seq_to_str(seq: Sequence) -> str:
    """Extract the sequence string from a Sequence object."""
    return seq.seqtext.as_str()


def make_run_codon_align(
    codon_align_fn: Callable[..., Any],
) -> Callable[..., tuple[str, str]]:
    """Factory: build a run_codon_align helper bound to a specific impl."""

    def run_codon_align(
        ref_str: str,
        seq_str: str,
        *,
        min_gap_distance: int = 30,
        window_size: int = 10,
        gap_placement_score: (dict[int, dict[tuple[int, int], int]] | None) = None,
        ref_start: int = 1,
        ref_end: int = -1,
    ) -> tuple[str, str]:
        if gap_placement_score is None:
            gap_placement_score = DEFAULT_GAP_PLACEMENT_SCORE

        refseq = make_sequence(ref_str, header='ref', seqid=0)
        seq = make_sequence(seq_str, header='seq', seqid=0)

        if ref_end <= 0:
            ref_end = refseq.seqtext.max_pos()

        ref_out, seq_out = codon_align_fn(
            refseq,
            seq,
            min_gap_distance,
            window_size,
            gap_placement_score,
            ref_start,
            ref_end,
        )
        return seq_to_str(ref_out), seq_to_str(seq_out)

    return run_codon_align
