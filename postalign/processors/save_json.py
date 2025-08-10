"""Processor generating NucAmino-style JSON reports."""

import orjson
# import math
import typer
from textwrap import indent
from typing import Tuple, List, Iterable, TypedDict, Optional, Dict
from more_itertools import chunked

from ..cli import cli
from ..utils import group_by_gene_codons, find_codon_trim_slice
from ..utils.codonutils import translate_codon
from ..models import RefSeqPair, NAPosition, Sequence, PositionFlag, Message
from ..processor import Processor, output_processor


GeneRanges = List[Tuple[str, List[Tuple[int, int]]]]


def gene_range_tuples_callback(
    ctx: typer.Context,
    param: typer.CallbackParam,
    value: Tuple[str, ...],
) -> List[Tuple[str, List[Tuple[int, int]]]]:
    one: str
    tuples: List[Tuple[str, List[Tuple[int, int]]]] = []
    cur_gene: str = ''
    cur_ranges: List[int] = []
    # use $$ as the end, it never goes into the results
    for one in value + ('$$', ):
        if one.isdigit():
            cur_ranges.append(int(one))
        else:
            if cur_gene and cur_ranges:
                if len(cur_ranges) % 2 != 0:
                    raise typer.BadParameter(
                        'Missing paired range values: <GENE>{}'
                        .format(cur_gene)
                    )
                cur_range_chunks: List[Tuple[int, int]] = []

                for refstart, refend in chunked(cur_ranges, 2):
                    if refstart < 1:
                        raise typer.BadParameter(
                            'argument <REF_START>:{} must be not less than 1'
                            .format(refstart)
                        )
                    if refend - 2 < refstart:
                        raise typer.BadParameter(
                            'no enough codon between arguments <REF_START>'
                            ':{} and <REF_END>:{}'
                            .format(refstart, refend)
                        )
                    cur_range_chunks.append((refstart, refend))
                tuples.append((cur_gene, cur_range_chunks))
            cur_gene = one
            cur_ranges = []
    return tuples


class CodonPair(TypedDict):
    Position: int
    RefCodonText: str
    CodonText: str
    RefAminoAcidText: str
    AminoAcidText: str
    InsertedCodonsText: str
    IsInsertion: bool
    IsDeletion: bool


class AlignedSite(TypedDict):
    PosAA: int
    PosNAs: List[Optional[int]]
    LengthNA: int


# XXX: remove total=False after upgraded to 3.11 (PEP 655)
class FrameShift(TypedDict, total=False):
    Position: int
    GapLength: int
    NucleicAcidsText: Optional[str]
    IsInsertion: int


class GeneReportInner(TypedDict, total=False):
    FirstAA: int
    LastAA: int
    AlignedSites: List[AlignedSite]
    Mutations: List[CodonPair]
    FrameShifts: List[FrameShift]


class GeneReport(TypedDict):
    Gene: str
    Report: GeneReportInner
    Error: str


class Payload(TypedDict):
    Name: str
    Modifiers: List[List[str]]
    GeneReports: List[GeneReport]
    Messages: List[Dict[str, str]]


@cli.command('save-json')
def save_json(
    gene_range_tuples: GeneRanges = typer.Argument(
        ..., callback=gene_range_tuples_callback,
    )
) -> Processor[Iterable[str]]:
    """Save as NucAmino style JSON reports

    <GENE_RANGE_TUPLES>:
        2n+1-tuples of
        <GENE> <REF_START1> <REF_END1> <REF_START2> <REF_END2> ...
        Use to calculate codon position within the protein/gene.
    """
    num_indent: int = 2

    @output_processor('save-json')
    def processor(
        iterator: Iterable[RefSeqPair],
        messages: List[Message]
    ) -> Iterable[str]:
        # TODO: MSA remap?
        refseq: Sequence
        seq: Sequence
        yield '[\n'
        for idx, (refseq, seq) in enumerate(iterator):
            reftext: List[NAPosition] = refseq.seqtext
            seqtext: List[NAPosition] = seq.seqtext
            seqmessages = [m for m in messages if m.seqid == seq.seqid]

            payload: Payload = {
                'Name': seq.description,
                'Modifiers': [[str(mod) for mod in mods]
                              for _, mods in seq.modifiers],
                'GeneReports': [],
                'Messages': [m.to_dict() for m in seqmessages]
            }

            gene_codons: List[
                Tuple[
                    str,
                    List[List[NAPosition]],
                    List[List[NAPosition]]
                ]
            ] = group_by_gene_codons(
                reftext, seqtext, gene_range_tuples)

            gene: str
            refcodons: List[List[NAPosition]]
            seqcodons: List[List[NAPosition]]
            for gene, refcodons, seqcodons in gene_codons:
                pos0: int
                refcd: List[NAPosition]
                seqcd: List[NAPosition]
                codonpairs: List[CodonPair] = []
                aligned_sites: List[AlignedSite] = []
                frameshifts: List[FrameShift] = []

                trim_slice: slice = find_codon_trim_slice(seqcodons)

                for pos0, (refcd, seqcd) in list(
                    enumerate(zip(refcodons, seqcodons))
                )[trim_slice]:
                    na: NAPosition
                    ins_fs_len: int
                    del_fs_len: int
                    if len(refcd) < 3:
                        # partial matched reference, skip
                        continue
                    nalen: int = NAPosition.count_nongaps(seqcd)
                    if 0 < nalen < 3:
                        ins_fs_len = 0
                        del_fs_len = 3 - nalen
                    else:
                        ins_fs_len = nalen % 3
                        del_fs_len = 0
                    codon_text: str = NAPosition.as_str(seqcd[:3])
                    if NAPosition.any_has_flag(
                        seqcd[:3], PositionFlag.TRIM_BY_SEQ
                    ):
                        continue
                    codonpairs.append({
                        'Position': pos0 + 1,
                        'RefCodonText': NAPosition.as_str(refcd[:3]),
                        'CodonText': codon_text,
                        'RefAminoAcidText': str(
                            translate_codon(refcd[:3]), 'ASCII'),
                        'AminoAcidText': str(
                            translate_codon(seqcd[:3]), 'ASCII'),
                        'InsertedCodonsText': NAPosition.as_str(
                            seqcd[3:len(seqcd) - ins_fs_len]
                        ),
                        'IsInsertion': len(refcd) > 5,
                        'IsDeletion': codon_text == '---'
                    })
                    posnas: List[Optional[int]] = []
                    for na in seqcd:
                        pos: int = na.pos
                        posnas.append(pos if pos > -1 else None)
                    aligned_sites.append({
                        'PosAA': pos0 + 1,  # reference location
                        'PosNAs': posnas,
                        'LengthNA': nalen
                    })
                    if ins_fs_len:
                        frameshifts.append({
                            'Position': pos0 + 1,
                            'GapLength': ins_fs_len,
                            'NucleicAcidsText':
                            NAPosition.as_str(seqcd[-ins_fs_len:]),
                            'IsInsertion': True
                        })
                    elif del_fs_len:
                        frameshifts.append({
                            'Position': pos0 + 1,
                            'GapLength': del_fs_len,
                            'IsInsertion': False
                        })
                mutations: List[CodonPair] = [
                    cd for cd in codonpairs
                    if cd['RefAminoAcidText'] != cd['AminoAcidText'] or
                    cd['IsInsertion']
                ]

                if codonpairs:
                    payload['GeneReports'].append({
                        'Gene': gene,
                        'Report': {
                            'FirstAA': codonpairs[0]['Position'],
                            'LastAA': codonpairs[-1]['Position'],
                            'AlignedSites': aligned_sites,
                            'Mutations': mutations,
                            'FrameShifts': frameshifts
                        },
                        'Error': ''
                    })
                else:
                    payload['GeneReports'].append({
                        'Gene': gene,
                        'Report': {},
                        'Error': 'Sequence is not aligned'
                    })

            text: bytes = orjson.dumps(payload, option=orjson.OPT_INDENT_2)
            if idx > 0:
                yield ',\n'
            yield indent(str(text, 'ASCII'), prefix=' ' * num_indent)
        yield '\n]\n'

    return processor
