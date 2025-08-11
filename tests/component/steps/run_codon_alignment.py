"""Execute codon alignment pipeline for component tests."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Callable, List
import inspect

import typer

if not hasattr(typer.Typer, "result_callback"):
    def result_callback(
        self: Any, *args: Any, **kwargs: Any
    ) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
            return func

        return decorator

    typer.Typer.result_callback = result_callback  # type: ignore[attr-defined]

if "multiple" not in inspect.signature(typer.Option).parameters:
    _orig_option = typer.Option

    def Option(*args: Any, **kwargs: Any) -> Any:  # type: ignore[override]
        kwargs.pop("multiple", None)
        return _orig_option(*args, **kwargs)

typer.Option = Option  # type: ignore[assignment]

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from postalign.parsers import minimap2  # noqa: E402
from postalign.processors.codon_alignment import codon_alignment  # noqa: E402
from postalign.models.sequence import NAPosition  # noqa: E402
from postalign.models import Message  # noqa: E402


def main() -> None:
    """Run codon alignment and write pairwise FASTA output.

    Expected arguments:
        seqs_path, ref_path, minimap2_bin, output_path, opts, start, end.
    """
    seqs_path = Path(sys.argv[1])
    ref_path = Path(sys.argv[2])
    minimap2_bin = sys.argv[3]
    output_path = Path(sys.argv[4])
    opts = sys.argv[5]
    start = int(sys.argv[6])
    end = int(sys.argv[7])

    messages: List[Message] = []
    with open(seqs_path) as seqs, open(ref_path) as ref:
        iterator = minimap2.load(
            seqs,
            ref,
            NAPosition,
            messages,
            minimap2_execute=[minimap2_bin, *opts.split()],
        )
        processor = codon_alignment(
            min_gap_distance=15, ref_start=start, ref_end=end
        )
        pairs = list(processor(iterator, messages))
    with open(output_path, "w") as out:
        for refseq, seq in pairs:
            out.write(f">{refseq.header}\n{refseq.seqtext_as_str}\n")
            out.write(f">{seq.header}\n{seq.seqtext_as_str}\n")


if __name__ == "__main__":
    main()
