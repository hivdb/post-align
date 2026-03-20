#! /usr/bin/env python
# -*- coding: UTF-8 -*-

# This file exists only for Cython extension builds.
# All metadata and config is in pyproject.toml.
# Run `make build-ext` (or `python setup.py build_ext --inplace`)
# to compile Cython extensions.

import sys

import setuptools

ext_modules = []

if 'build_ext' in sys.argv:
    from Cython.Build import cythonize  # type: ignore
    from setuptools.extension import Extension

    extensions = [
        # --- shared utilities ---
        Extension(
            name='postalign.utils.cigar',
            sources=['postalign/utils/cigar.py']
        ),
        Extension(
            name='postalign.models.na_position',
            sources=['postalign/models/na_position.py']
        ),
        Extension(
            name='postalign.models._sequence',
            sources=['postalign/models/_sequence.py']
        ),
        Extension(
            name='postalign.utils.group_by_codons',
            sources=['postalign/utils/group_by_codons.py']
        ),
        Extension(
            name='postalign.utils.codonutils',
            sources=['postalign/utils/codonutils.py']
        ),
        Extension(
            name='postalign.utils.blosum62',
            sources=['postalign/utils/blosum62.py']
        ),
        Extension(
            name='postalign.utils.iupac',
            sources=['postalign/utils/iupac.py']
        ),
        # --- T1: original ---
        Extension(
            name='postalign.processors.codon_alignment',
            sources=['postalign/processors/codon_alignment.py']
        ),
        # --- T2: optimized Python + D ---
        Extension(
            name='postalign.processors.codon_alignment_optimized',
            sources=['postalign/processors/codon_alignment_optimized.py']
        ),
        # --- T3: Cython byte-array + D ---
        Extension(
            name='postalign.processors.codon_alignment_cython',
            sources=['postalign/processors/codon_alignment_cython.py']
        ),
        # --- T4: Rust wrapper (thin Python layer) ---
        Extension(
            name='postalign.processors.codon_alignment_rust',
            sources=['postalign/processors/codon_alignment_rust.py']
        ),
    ]
    ext_modules = cythonize(
        extensions,
        compiler_directives={
            'language_level': '3',
            'profile': False,
            'linetrace': False,
            'cdivision': True,
        }
    )

if __name__ == '__main__':
    setuptools.setup(ext_modules=ext_modules)
