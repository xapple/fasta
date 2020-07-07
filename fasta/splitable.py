#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, sys, math, shutil

# Internal modules #
from fasta import FASTA

# Third party modules #

###############################################################################
class SplitableFASTA(FASTA):
    """
    A FASTA file which you can split into chunks. Either you give the number
    of parts you want to generate, or you can give a target size in bytes for
    each part.
    """

    def __init__(self, path, num_parts=None, part_size=None, base_dir=None):
        # Basic #
        self.path = path
        # Directory #
        if base_dir is None: self.base_dir = path + '.parts/'
        else:                self.base_dir = base_dir
        # Num parts #
        if num_parts is not None: self.num_parts = num_parts
        # Special module #
        import humanfriendly
        # Evaluate size #
        if part_size is not None:
            self.bytes_target = humanfriendly.parse_size(part_size)
            self.num_parts = int(math.ceil(self.count_bytes / self.bytes_target))
        # Make parts #
        self.make_name = lambda i: self.base_dir + "%03d/part.fasta" % i
        self.parts = [FASTA(self.make_name(i)) for i in range(1, self.num_parts+1)]
        # Give a number to each part #
        for i, part in enumerate(self.parts): part.num = i

    @property
    def status(self):
        """Has the splitting been done already?"""
        if all(os.path.exists(p.path) for p in self.parts): return True
        return False

    def run(self):
        # Clean up #
        for i in range(1, sys.maxsize):
            dir_path = self.base_dir + "%03d/" % i
            if os.path.exists(dir_path): shutil.rmtree(dir_path)
            else: break
        # Case only one part #
        if len(self.parts) == 1:
            self.parts[0].directory.create(safe=True)
            self.link_to(self.parts[0])
            return
        # Compute number of sequences #
        self.seqs_per_part = int(math.floor(self.count / self.num_parts))
        # Prepare #
        for part in self.parts: part.create()
        # Do the job #
        seqs = self.parse()
        for part in self.parts:
            for i in range(self.seqs_per_part):
                part.add_seq(seqs.next())
        # The final sequences go to the last part #
        for seq in seqs: part.add_seq(seq)
        # Clean up #
        for part in self.parts: part.close()