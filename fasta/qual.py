#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Internal modules #
from fasta import FASTA

###############################################################################
class QualFile(FASTA):
    """A single QUAL file somewhere in the filesystem."""

    format = 'qual'
    ext    = 'qual'