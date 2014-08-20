# Futures #
from __future__ import division

# Built-in modules #
import re

# Internal modules #
from fasta import FASTA
from plumbing.cache import property_cached

###############################################################################
class SizesFASTA(FASTA):
    """A FASTA file with size annotation affecting the count"""

    @property_cached
    def count(self):
        get_size = lambda x: int(re.findall("size=([0-9]+)", x)[0])
        sizes = (get_size(r.description) for r in self)
        return sum(sizes)