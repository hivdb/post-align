#! /usr/bin/env python

from pathlib import Path

import setuptools
from Cython.Build import cythonize  # type: ignore
from setuptools.extension import Extension

extensions = [
    Extension(
        name="postalign.utils.cigar",
        sources=["postalign/utils/cigar.py"],
    ),
    Extension(
        name="postalign.models.na_position",
        sources=["postalign/models/na_position.py"],
    ),
    Extension(
        name="postalign.models._sequence",
        sources=["postalign/models/_sequence.py"],
    ),
    Extension(
        name="postalign.processors.codon_alignment",
        sources=["postalign/processors/codon_alignment.py"],
    ),
    Extension(
        name="postalign.utils.group_by_codons",
        sources=["postalign/utils/group_by_codons.py"],
    ),
    Extension(
        name="postalign.utils.codonutils",
        sources=["postalign/utils/codonutils.py"],
    ),
    Extension(
        name="postalign.utils.blosum62",
        sources=["postalign/utils/blosum62.py"],
    ),
    Extension(
        name="postalign.utils.iupac",
        sources=["postalign/utils/iupac.py"],
    ),
]


def strip_comments(line: str) -> str:
    """Remove comments and index options from a requirement line."""

    if line.startswith("-i "):
        return ""
    return line.split("#", 1)[0].strip()


def req(filename: str) -> list[str]:
    """Load dependencies from *filename* and return a list of requirements."""

    with Path(filename).open() as fp:
        requires: set[str] = {strip_comments(ln) for ln in fp.readlines()}
        requires.discard("")
    return list(requires)


if __name__ == "__main__":
    setuptools.setup(
        ext_modules=cythonize(
            extensions,
            compiler_directives={
                "language_level": "3",
                "profile": False,
                "linetrace": False,
            },
        ),  # type: ignore[no-untyped-call]
        install_requires=req("requirements.txt"),
    )
