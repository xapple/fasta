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
import inspect

# First party modules #
from autopaths import Path

# Third party modules #

# Constants #
this_file = Path((inspect.stack()[0])[1])
this_dir  = this_file.directory

###############################################################################
# --------------------------------- Writing --------------------------------- #
def test_write(tmp_path):
    # Import #
    from fasta import FASTQ
    # Create object #
    fastq = this_dir + "data/seqs.fastq"
    fastq = FASTQ(fastq)
    # Create output #
    out = str(tmp_path / 'new.fastq')
    out = FASTQ(out)
    # Copy sequences #
    out.write(iter(fastq))
    # Check length #
    assert len(out) == 1401

def test_write_gz(tmp_path):
    # Import #
    from fasta import FASTQ
    # Create object #
    fastq = this_dir + "data/seqs.fastq.gz"
    fastq = FASTQ(fastq)
    # Create output #
    out = str(tmp_path / 'new.fastq.gz')
    out = FASTQ(out)
    # Copy sequences #
    out.write(iter(fastq))
    # Check length #
    assert len(out) == 1401
