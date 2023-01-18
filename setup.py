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
from platform import system

# Load the contents of the README file #
this_dir = path.abspath(path.dirname(__file__))
readme_path = path.join(this_dir, 'README.md')
with open(readme_path, encoding='utf-8') as handle: readme = handle.read()

install_requires = ['autopaths>=1.5.0', 'plumbing>=2.10.4', 'biopython', 'tqdm']
if system() == 'Windows': 
    install_requires.append('pbs3')
else:
    install_requires.append('sh')

# Call setup #
setup(
    name             = 'fasta',
    version          = '2.2.14',
    description      = 'The fasta python package enables you to deal with '
                       'biological sequence files easily.',
    license          = 'MIT',
    url              = 'http://github.com/xapple/fasta/',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
    packages         = find_namespace_packages(),
    install_requires = install_requires,
    extras_require   = {'graphs':  ['numpy', 'matplotlib'],
                        'primers': ['regex']},
    python_requires  = ">=3.8",
    long_description = readme,
    long_description_content_type = 'text/markdown',
    include_package_data = True,
)
