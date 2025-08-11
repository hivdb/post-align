#! /usr/bin/env python
"""Build configuration for Cython extensions."""

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
    )
