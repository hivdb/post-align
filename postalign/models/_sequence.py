import cython  # type: ignore

from .aa_position import AAPosition
from .na_position import NAPosition, NAPositionList

type Position = NAPosition | AAPosition


class _SkipValidationSentinel:
    """Sentinel type used to bypass sequence validation."""


SKIP_VALIDATION = _SkipValidationSentinel()

VALID_NOTATIONS: dict[type[NAPosition] | type[AAPosition], set[int]] = {NAPosition: set(b'ACGTUWSMKRYBDHVN.-')}


@cython.ccall
def sanitize_sequence(
    seqtext: list[NAPosition] | list[AAPosition] | NAPositionList,
    seqtype: type[NAPosition] | type[AAPosition],
    header: str,
    skip_invalid: bool | _SkipValidationSentinel,
) -> list[NAPosition] | list[AAPosition] | NAPositionList:
    if seqtype == AAPosition:
        raise NotImplementedError('Amino acid is not yet supported')

    if isinstance(seqtext, NAPositionList):
        # NAPositionList is inherently all-NAPosition; just validate notations
        pass
    elif not all(isinstance(one, NAPosition) for one in seqtext):
        raise ValueError("seqtext must be a list of NAPosition instances when seqtype is 'NAPosition'")

    valid_notations: set[int] = VALID_NOTATIONS[seqtype]
    valids: list[NAPosition] = []
    invalids: set[int] = set()
    for one in seqtext:
        if isinstance(one, NAPosition):
            if one.notation in valid_notations:
                valids.append(one)
            else:
                invalids.add(one.notation)
    if invalids and skip_invalid:
        seqtext = NAPositionList.from_list(valids) if isinstance(seqtext, NAPositionList) else valids
    elif invalids:
        raise ValueError(
            'sequence {} contains invalid notation(s) ({})while skip_invalid=False'.format(
                header, str(bytes(sorted(invalids)), 'ASCII')
            )
        )
    return seqtext
