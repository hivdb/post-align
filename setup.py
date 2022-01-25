#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import re
import ast
import setuptools
from Cython.Build import cythonize  # type: ignore
from typing import Optional, List, Set
from setuptools.extension import Extension

version: str
_version_re: re.Pattern = re.compile(r'VERSION\s+=\s+(.*)')

extensions = [
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
        name='postalign.processors.codon_alignment',
        sources=[
            'postalign/processors/codon_alignment.py'
        ]
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
    )
]

with open('postalign/version.py', 'rb') as f:
    match: Optional[re.Match] = _version_re.search(
        f.read().decode('utf-8'))
    if match:
        version = str(ast.literal_eval(match.group(1)))
    else:
        raise ImportError('Unable to import version from postalign.version')


def strip_comments(line: str) -> str:
    if line.startswith('-i '):
        return ''
    return line.split('#', 1)[0].strip()


def req(filename: str) -> List[str]:
    with open(os.path.join(os.getcwd(), filename)) as fp:
        requires: Set[str] = set([strip_comments(ln) for ln in fp.readlines()])
        requires -= set([''])
    return list(requires)


if __name__ == '__main__':
    setuptools.setup(
        name='post-align',
        version=version,
        url="https://github.com/hivdb/post-align",
        author='Philip Tzou',
        author_email="philiptz@stanford.edu",
        description=(
            'A post alignment toolkit for refining '
            'pairwise/multiple alignment sequences.'
        ),
        packages=[
            'postalign', 'postalign/models',
            'postalign/parsers', 'postalign/processors',
            'postalign/utils'
        ],
        package_data={"postalign": ["py.typed"]},
        install_requires=req('requirements.txt'),
        ext_modules=cythonize(
            extensions,
            compiler_directives={
                'language_level': '3',
                'profile': False,
                'linetrace': False
            }
        ),
        # tests_require=reqs('test-requirements.txt'),
        # include_package_data=True,
        entry_points={'console_scripts': [
            'postalign = postalign.entry:cli',
        ]},
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Topic :: Scientific/Engineering :: Bio-Informatics'],
        # test_suite="nose.collector",
        zip_safe=True)
