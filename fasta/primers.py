#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import re

# Internal modules #
from plumbing.common import GenWithLength
from plumbing.color  import Color

# Third party modules #
from Bio.Seq import Seq

# Constants #
iupac = {'A':'A',    'G':'G',   'T':'T',   'C':'C',
         'M':'AC',   'R':'AG',  'W':'AT',  'S':'CG',   'Y':'CT',   'K':'GT',
         'V':'ACG',  'H':'ACT', 'D':'AGT', 'B':'CGT',
         'X':'ACGT', 'N':'ACGT'}

# Function #
iupac_pattern = lambda seq: ''.join(['[' + iupac[char] + ']' for char in seq])

################################################################################
def parse_primers(fasta, primers=None, mismatches=0, revcompl=False):
    """
    This functions takes a FASTA instance as first parameter.
    The primers will be loaded from that FASTA if not specified.
    """
    # Default primers #
    if primers is None: primers = fasta.primers
    # Special module #
    import regex
    # Case straight #
    if not revcompl:
        fwd_regex = regex.compile("(%s){s<=%i}" % (primers.fwd_pattern, mismatches))
        rev_regex = regex.compile("(%s){s<=%i}" % (primers.rev_pattern, mismatches))
        generator = (ReadWithPrimers(r, fwd_regex, rev_regex) for r in fasta.parse())
    # Case revcompl #
    if revcompl:
        fwd_regex = regex.compile("(%s){s<=%i}" % (primers.fwd_pattern,          mismatches))
        rev_regex = regex.compile("(%s){s<=%i}" % (primers.rev_pattern_revcompl, mismatches))
        generator = (ReadWithPrimersRevCompl(r, fwd_regex, rev_regex) for r in fasta.parse())
    # Return #
    return GenWithLength(generator, len(fasta))

###############################################################################
class TwoPrimers:
    """A container for the two primers of a sample."""

    def __len__(self): return 2

    def __init__(self, fwd_str, rev_str):
        # Strings #
        self.fwd_str = fwd_str
        self.rev_str = rev_str
        # Lengths #
        self.fwd_len = len(self.fwd_str)
        self.rev_len = len(self.rev_str)
        # Sequences #
        self.fwd_seq = Seq(self.fwd_str)
        self.rev_seq = Seq(self.rev_str)
        # Search patterns #
        self.fwd_pattern = iupac_pattern(self.fwd_seq)
        self.rev_pattern = iupac_pattern(self.rev_seq) # Don't add reverse complement here, use option instead
        # Search patterns reverse complemented #
        self.fwd_pattern_revcompl = iupac_pattern(self.fwd_seq.reverse_complement())
        self.rev_pattern_revcompl = iupac_pattern(self.rev_seq.reverse_complement())
        # Search expression without mismatches #
        self.fwd_regex = re.compile(self.fwd_pattern)
        self.rev_regex = re.compile(self.rev_pattern)
        # Uracil instead of thymine #
        self.fwd_regex_uracil = re.compile(self.fwd_pattern.replace('T', 'U'))
        self.rev_regex_uracil = re.compile(self.rev_pattern.replace('T', 'U'))

###############################################################################
class ReadWithPrimers:
    def __init__(self, read, fwd_regex, rev_regex):
        self.read          = read
        self.fwd_match     = fwd_regex.search(str(read.seq))
        self.rev_match     = rev_regex.search(str(read.seq))
        self.fwd_start_pos = self.fwd_match.start() if self.fwd_match else None
        self.rev_start_pos = self.rev_match.start() if self.rev_match else None
        self.fwd_end_pos   = self.fwd_match.end()   if self.fwd_match else None
        self.rev_end_pos   = self.rev_match.end()   if self.rev_match else None

    @property
    def pretty(self):
        # The string #
        seq  = self.read.seq._data
        fwds = self.fwd_start_pos
        fwde = self.fwd_end_pos
        revs = self.rev_start_pos
        reve = self.rev_end_pos
        # No matches #
        if not self.fwd_match and not self.rev_match:
            return seq + '\n'
        # One matches #
        if self.fwd_match and not self.rev_match:
            return seq[0:fwds]    + Color.red + \
                   seq[fwds:fwde] + Color.end + \
                   seq[fwde:]     + '\n'
        # One matches #
        if not self.fwd_match and self.rev_match:
            return seq[0:revs]    + Color.red + \
                   seq[revs:reve] + Color.end + \
                   seq[reve:]     + '\n'
        # Two matches #
        if self.fwd_match and self.rev_match:
            return seq[0:fwds]    + Color.red + \
                   seq[fwds:fwde] + Color.end + \
                   seq[fwde:revs] + Color.red + \
                   seq[revs:reve] + Color.end + \
                   seq[reve:]     + '\n'

###############################################################################
class ReadWithPrimersRevCompl:
    """Here the reverse primer start and end values will be negative."""

    def __init__(self, read, fwd_regex, rev_regex):
        self.read          = read
        self.fwd_match     = fwd_regex.search(str(read.seq))
        self.rev_match     = rev_regex.search(str(read.seq))
        self.fwd_start_pos = self.fwd_match.start()             if self.fwd_match else None
        self.rev_start_pos = self.rev_match.end() - len(read)   if self.rev_match else None
        self.fwd_end_pos   = self.fwd_match.end()               if self.fwd_match else None
        self.rev_end_pos   = self.rev_match.start() - len(read) if self.rev_match else None