#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

This file contains tests to be run automatically with the pytest
executable. To run all the tests just type the following on your
terminal from the repository root:

    $ python3 -m pip install --upgrade --user pytest
    $ pytest
"""

# Built-in modules #
import sys, inspect

# First party modules #
from autopaths import Path

# Third party modules #

# Constants #
this_file = Path((inspect.stack()[0])[1])
this_dir  = this_file.directory

###############################################################################
# -------------------------------- Counting --------------------------------- #
def test_count():
    valid_file = this_dir + "data/seqs.fastq"
    from fasta import FASTQ
    fastq = FASTQ(valid_file)
    assert len(fastq) == 1401

def test_count_gz():
    valid_file = this_dir + "data/seqs.fastq.gz"
    from fasta import FASTQ
    fastq = FASTQ(valid_file)
    assert len(fastq) == 1401
