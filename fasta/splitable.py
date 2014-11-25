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

    def __init__(self, path, part_size, base_dir=None):
        # Basic #
        self.path = path
        # Directory #
        if base_dir is None: self.base_dir = path + '.parts/'
        else:                self.base_dir = base_dir
        # Evaluate size #
        self.bytes_target = humanfriendly.parse_size(part_size)
        # Chose number of parts #
        self.num_parts = int(math.ceil(self.count_bytes / self.bytes_target))
        # Make parts #
        self.make_name = lambda i: self.base_dir + "%03d/part.fasta" % i
        self.parts = [FASTA(self.make_name(i)) for i in range(self.num_parts)]

    @property
    def status(self):
        if all(os.path.exists(p.path) for p in self.parts): return 'splitted'
        return False

    def split(self):
        # Clean #
        for i in xrange(sys.maxint):
            if os.path.exists(self.make_name(i)): os.remove(self.make_name(i))
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
            for i in xrange(self.seqs_per_part): part.add_read(seqs.next())
        for seq in seqs: part.add_read(seq)
        # Clean up #
        for part in self.parts: part.close()