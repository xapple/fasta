# Futures #
from __future__ import division

# Built-in modules #
import os, sys, math

# Internal modules #
from fasta import FASTA

# Third party modules #
import humanfriendly

###############################################################################
class SplitableFASTA(FASTA):
    """A FASTA file which you can split into chunks. Either you give the number
    of parts you want to generate, or you can give a target size in bytes for
    each part."""

    def __init__(self, path, num_parts=None, part_size=None, base_dir=None):
        # Basic #
        self.path = path
        # Directory #
        if base_dir is None: self.base_dir = path + '.parts/'
        else: self.base_dir = base_dir
        if not os.path.exists(self.base_dir): os.makedirs(self.base_dir)
        # Num parts #
        if num_parts is not None: self.num_parts = num_parts
        # Evaluate size #
        if part_size is not None:
            self.bytes_target = humanfriendly.parse_size(part_size)
            self.num_parts = int(math.ceil(self.count_bytes / self.bytes_target))
        # Make parts #
        self.parts = [FASTA(self.make_part_paths(i)) for i in range(self.num_parts)]
        # Give numbers #
        for i, part in enumerate(self.parts): part.num = i

    def make_part_paths(self, i):
        """Generate the paths for the different parts"""
        return self.base_dir + "%s%03d.fasta" % (self.prefix, i)

    @property
    def status(self):
        """Has the splitting been done already ?"""
        if all(os.path.exists(p.path) for p in self.parts): return 'splitted'
        return False

    def split(self):
        # Clean #
        for i in xrange(sys.maxint):
            if os.path.exists(self.make_part_paths(i)): os.remove(self.make_part_paths(i))
            else: break
        # Case only one part #
        if len(self.parts) == 1:
            os.symlink(self.path, self.parts[0].path)
            return
        # Compute number of sequences #
        self.seqs_per_part = int(math.floor(self.count / self.num_parts))
        # Prepare #
        for part in self.parts: part.create()
        # Do the job #
        seqs = self.parse()
        for part in self.parts:
            for i in xrange(self.seqs_per_part): part.add_seq(seqs.next())
        for seq in seqs: part.add_seq(seq)
        # Clean up #
        for part in self.parts: part.close()