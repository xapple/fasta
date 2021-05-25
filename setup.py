#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Imports #
from setuptools import setup, find_namespace_packages
from os import path

# Load the contents of the README file #
this_dir = path.abspath(path.dirname(__file__))
readme_path = path.join(this_dir, 'README.md')
with open(readme_path, encoding='utf-8') as handle: readme = handle.read()

# Call setup #
setup(
    name             = 'fasta',
    version          = '2.2.11',
    description      = 'The fasta python package enables you to deal with '
                       'biological sequence files easily.',
    license          = 'MIT',
    url              = 'http://github.com/xapple/fasta/',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
    packages         = find_namespace_packages(),
    install_requires = ['plumbing>=2.10.4', 'autopaths>=1.4.12', 'biopython',
                        'sh', 'tqdm'],
    extras_require   = {'graphs':  ['numpy', 'matplotlib'],
                        'primers': ['regex']},
    python_requires  = ">=3.8",
    long_description = readme,
    long_description_content_type = 'text/markdown',
    include_package_data = True,
)
