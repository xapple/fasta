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
from plumbing.color  import Color

# Third party modules #
import Bio, regex
from Bio.Seq import Seq

# Constants #
iupac = {'A':'A',    'G':'G',   'T':'T',   'C':'C',
         'M':'AC',   'R':'AG',  'W':'AT',  'S':'CG',   'Y':'CT',   'K':'GT',
         'V':'ACG',  'H':'ACT', 'D':'AGT', 'B':'CGT',
         'X':'ACGT', 'N':'ACGT'}

# Function to create a regex pattern from a sequence #
iupac_pattern = lambda seq: ''.join(['[' + iupac[char] + ']' for char in seq])

###############################################################################
class TwoPrimers:
    """
    A container for the two primers of a sample.
    Has methods for generating regexes to search for these primers.
    """

    def __len__(self): return 2

    def __init__(self, fwd_str, rev_str):
        # Original strings #
        self.fwd_str = fwd_str
        self.rev_str = rev_str
        # Lengths in base pairs #
        self.fwd_len = len(self.fwd_str)
        self.rev_len = len(self.rev_str)
        # Sequences as biopython objects #
        self.fwd_seq = Bio.Seq.Seq(self.fwd_str)
        self.rev_seq = Bio.Seq.Seq(self.rev_str)
        # Create search patterns in regex syntax #
        self.fwd_pat = iupac_pattern(self.fwd_seq)
        self.rev_pat = iupac_pattern(self.rev_seq)
        # Reverse complemented sequences #
        self.fwd_revcomp = self.fwd_seq.reverse_complement()
        self.rev_revcomp = self.rev_seq.reverse_complement()
        # Search patterns when reverse complemented #
        self.fwd_pat_revcomp = iupac_pattern(self.fwd_revcomp)
        self.rev_pat_revcomp = iupac_pattern(self.rev_revcomp)
        # Simple search expression (without any mismatches authorized yet) #
        self.fwd_search = re.compile(self.fwd_pat)
        self.rev_search = re.compile(self.rev_pat)

    def make_regex(self, pat, mismatches):
        """Complex search expression with mismatches this time."""
        return regex.compile("(%s){s<=%i}" % (pat, mismatches))

    def make_fwd_regex(self, mismatches):
        return self.make_regex(self.fwd_pat, mismatches)

    def make_rev_regex(self, mismatches):
        return self.make_regex(self.rev_pat, mismatches)

    def make_fwd_revcompl_regex(self, mismatches):
        return self.make_regex(self.fwd_pat_revcomp, mismatches)

    def make_rev_revcompl_regex(self, mismatches):
        return self.make_regex(self.rev_pat_revcomp, mismatches)

###############################################################################
class PrimersRegexes:
    """
    A container for the regular expression search patterns
    that enable us to find primers inside a sequence.
    These regexes depend on the number of mismatches authorized.
    """

    def __init__(self, primers, mismatches):
        """
        We need to know the primers and the number of mismatches tolerated
        in the search.
        """
        # Base attributes #
        self.primers    = primers
        self.mismatches = mismatches
        # Search patterns #
        self.fwd    = primers.make_fwd_regex(mismatches)
        self.rev    = primers.make_rev_regex(mismatches)
        # Search patterns reverse complemented #
        self.fwd_rc = primers.make_fwd_revcompl_regex(mismatches)
        self.rev_rc = primers.make_rev_revcompl_regex(mismatches)

###############################################################################
class ReadWithPrimers:
    def __init__(self, read, regexes):
        """
        Uses regex patterns to search the given read.
        Records the start and end positions of primers if they are found.
        Both the forward and reverse primers are searched for.
        Both the original sequences and their reverse complements are
        searched for, in case the read is in the opposite direction.
        """
        # The read itself #
        self.read = read
        # The sequence as a string #
        self.seq = str(read.seq)
        # Searches #
        self.fwd    = regexes.fwd.search(self.seq)
        self.rev    = regexes.rev.search(self.seq)
        self.fwd_rc = regexes.fwd_rc.search(self.seq)
        self.rev_rc = regexes.rev_rc.search(self.seq)
        # Positions found in standard search #
        self.fwd_srt = self.fwd.start() if self.fwd else None
        self.fwd_end = self.fwd.end()   if self.fwd else None
        self.rev_srt = self.rev.start() if self.rev else None
        self.rev_end = self.rev.end()   if self.rev else None
        # Positions found in reverse complement search #
        self.fwd_rc_srt = self.fwd_rc.start() if self.fwd_rc else None
        self.fwd_rc_end = self.fwd_rc.end()   if self.fwd_rc else None
        self.rev_rc_srt = self.rev_rc.start() if self.rev_rc else None
        self.rev_rc_end = self.rev_rc.end()   if self.rev_rc else None

    @property
    def pretty_visualization(self):
        """
        This property is useful for debugging.
        It will return a nicely formatted string showing the original read
        with all primers found highlighted with bash color codes.
        """
        # Make a copy of the read for convenience #
        seq = self.seq
        # Initialize output #
        out = ""
        # Iterate over every position in the original sequence #
        for i, nuc in enumerate(seq):
            if i == self.fwd_srt:    out += Color.b_grn
            if i == self.rev_srt:    out += Color.grn
            if i == self.fwd_rc_srt: out += Color.red
            if i == self.rev_rc_srt: out += Color.b_red
            if i == self.fwd_end:    out += Color.end
            if i == self.rev_end:    out += Color.end
            if i == self.fwd_rc_end: out += Color.end
            if i == self.rev_rc_end: out += Color.end
            out += nuc
        # Summary of found positions #
        summary = f"""
        Forward start:           {self.fwd_srt}
        Forward end:             {self.fwd_end}
        Reverse start:           {self.rev_srt}
        Reverse end:             {self.rev_end}
        Forward revcompl start:  {self.fwd_rc_srt}
        Forward revcompl end:    {self.fwd_rc_end}
        Reverse revcompl start:  {self.rev_rc_srt}
        Reverse revcompl end:    {self.rev_rc_end}
        """
        # Return #
        return summary + out + '\n'
