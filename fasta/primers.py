# Built-in modules #
import re

# Internal modules #
from plumbing.common import GenWithLength

# Third party modules #
import regex
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

# Constants #
iupac = {'A':'A', 'G':'G', 'T':'T', 'C':'C', 'M':'AC', 'R':'AG', 'W':'AT', 'S':'CG', 'Y':'CT', 'K':'GT', 'V':'ACG', 'H':'ACT', 'D':'AGT', 'B':'CGT', 'X':'ACGT', 'N':'ACGT'}

################################################################################
def parse_primers(self, primers=None, mismatches=0):
    # Default primers #
    if primers is None: primers = self.primers
    # Case no mismatches #
    if mismatches == 0:
        generator = (ReadWithPrimers(r, primers) for r in self.parse())
        return GenWithLength(generator, len(self))
    # Case with mismatches #
    else:
        fwd_regex = regex.compile("(%s){s<=%i}" % (primers.fwd_pattern, mismatches))
        rev_regex = regex.compile("(%s){s<=%i}" % (primers.rev_pattern, mismatches))
        generator = (ReadWithPrimersMissmatch(r, primers, fwd_regex, rev_regex) for r in self.parse())
        return GenWithLength(generator, len(self))

###############################################################################
class TwoPrimers(object):
    """A container for the two primers of a sample"""

    def __repr__(self): return '<%s object for pool %s>' % (self.__class__.__name__, self.parent.id_name)
    def __len__(self): return 2

    def __init__(self, fwd_str, rev_str, fwd_name=None, rev_name=None):
        # Names #
        self.fwd_name = fwd_name
        self.rev_name = rev_name
        # Strings #
        self.fwd_str = fwd_str
        self.rev_str = rev_str
        # Lengths #
        self.fwd_len = len(self.fwd_str)
        self.rev_len = len(self.rev_str)
        # Sequences #
        self.fwd_seq = Seq(self.fwd_str, IUPAC.ambiguous_dna)
        self.rev_seq = Seq(self.rev_str, IUPAC.ambiguous_dna)
        # Search patterns #
        self.fwd_pattern = ''.join(['[' + iupac[char] + ']' for char in self.fwd_seq])
        self.rev_pattern = ''.join(['[' + iupac[char] + ']' for char in self.rev_seq.reverse_complement()])
        # Search expression #
        self.fwd_regex = re.compile(self.fwd_pattern)
        self.rev_regex = re.compile(self.rev_pattern)
        # Uracil instead of thymine #
        self.fwd_regex_uracil = re.compile(self.fwd_pattern.replace('T', 'U'))
        self.rev_regex_uracil = re.compile(self.rev_pattern.replace('T', 'U'))

###############################################################################
class ReadWithPrimers(object):
    def __init__(self, read, primers):
        self.read = read
        self.fwd_match     = primers.fwd_regex.search(str(read.seq))
        self.rev_match     = primers.rev_regex.search(str(read.seq))
        self.fwd_start_pos = self.fwd_match.start() if self.fwd_match else None
        self.rev_start_pos = self.rev_match.end() - len(read) if self.rev_match else None
        self.fwd_end_pos   = self.fwd_match.end() if self.fwd_match else None
        self.rev_end_pos   = self.rev_match.start() - len(read) if self.rev_match else None

###############################################################################
class ReadWithPrimersMissmatch(object):
    def __init__(self, read, primers, fwd_regex, rev_regex):
        self.read = read
        self.fwd_match     = fwd_regex.search(str(read.seq))
        self.rev_match     = rev_regex.search(str(read.seq))
        self.fwd_start_pos = self.fwd_match.start() if self.fwd_match else None
        self.rev_start_pos = self.rev_match.end() - len(read) if self.rev_match else None
        self.fwd_end_pos   = self.fwd_match.end() if self.fwd_match else None
        self.rev_end_pos   = self.rev_match.start() - len(read) if self.rev_match else None

