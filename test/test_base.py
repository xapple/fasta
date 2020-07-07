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
import pytest

# Constants #
this_file = Path((inspect.stack()[0])[1])
this_dir  = this_file.directory

###############################################################################
# --------------------------------- Import ---------------------------------- #
def test_package_import():
    import fasta
    return fasta

# ------------------------------ Validation --------------------------------- #
@pytest.mark.skipif(sys.platform != 'linux', reason="Can only run on Linux.")
class TestValidation:

    def test_validation():
        valid_file = this_dir + "data/seqs.fastq"
        from fasta import FASTQ
        fastq = FASTQ(valid_file)
        assert fastq.validator()

    def test_invalid():
        invalid_file = this_dir + "data/invalid.fastq"
        from fasta import FASTQ
        fastq = FASTQ(invalid_file)
        assert not fastq.validator()
