#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import re
import ast
import setuptools

_version_re = re.compile(r'VERSION\s+=\s+(.*)')

with open('postalign/version.py', 'rb') as f:
    version = str(ast.literal_eval(_version_re.search(
        f.read().decode('utf-8')).group(1)))


def strip_comments(line):
    if line.startswith('-i '):
        return ''
    return line.split('#', 1)[0].strip()


def req(filename):
    with open(os.path.join(os.getcwd(), filename)) as fp:
        requires = set([strip_comments(ln) for ln in fp.readlines()])
        requires -= set([''])
    return list(requires)


setup_params = dict(
    name="post-align",
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
    install_requires=req('requirements.txt'),
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

if __name__ == '__main__':
    setuptools.setup(**setup_params)
