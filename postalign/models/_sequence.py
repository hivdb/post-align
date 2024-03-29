import cython  # type: ignore
from typing import Union, List, Dict, Type, Set

from .na_position import NAPosition
from .aa_position import AAPosition

Position = Union[NAPosition, AAPosition]

SKIP_VALIDATION = object()

VALID_NOTATIONS: Dict[Type[Position], Set[int]] = {
    NAPosition: set(b'ACGTUWSMKRYBDHVN.-')
}


@cython.ccall
@cython.returns(list)
def sanitize_sequence(
    seqtext: List[Position],
    seqtype: Type[Position],
    header: str,
    skip_invalid: Union[bool, object]
) -> List[Position]:
    if (
        seqtype == AAPosition or
        any([isinstance(one, AAPosition) for one in seqtext])
    ):
        raise NotImplementedError('Amino acid is not yet supported')

    if (
        seqtype != NAPosition and
        not all([isinstance(one, NAPosition) for one in seqtext])
    ):
        raise ValueError(
            "seqtext must be a list of NAPosition instances "
            "when seqtype is 'NAPosition'"
        )

    valid_notations: set[int] = VALID_NOTATIONS[seqtype]
    valids: List[Position] = []
    invalids: Set[int] = set()
    for one in seqtext:
        if one.notation in valid_notations:
            valids.append(one)
        else:
            invalids.add(one.notation)
    if invalids and skip_invalid:
        seqtext = valids
    elif invalids:
        raise ValueError(
            'sequence {} contains invalid notation(s) ({})'
            'while skip_invalid=False'
            .format(
                header,
                str(bytes(sorted(invalids)), 'ASCII')
            )
        )
    return seqtext
