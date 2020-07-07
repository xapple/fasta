#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
from six.moves import zip as izip

# Internal modules #
from fasta import FASTA, FASTQ
from plumbing.common import isubsample, GenWithLength
from plumbing.cache import property_cached

# Third party modules #
from tqdm import tqdm

###############################################################################
class PairedFASTA:
    """Read and write FASTA file pairs without using too much RAM."""

    format = 'fasta'

    def __len__(self):     return self.count
    def __iter__(self):    return self.parse()
    def __bool__(self): return bool(self.fwd) and bool(self.rev)
    def __repr__(self): return '<%s object on "%s" and "%s">' % \
                        (self.__class__.__name__, self.fwd.path, self.rev.path)

    def __enter__(self): return self.create()
    def __exit__(self, exc_type, exc_value, traceback): self.close()

    @property
    def exists(self): return self.fwd.exists and self.rev.exists

    def __init__(self, fwd, rev, parent=None):
        # FASTA objects #
        self.fwd = FASTA(fwd)
        self.rev = FASTA(rev)
        # Extra #
        self.gzipped = self.fwd.gzipped
        self.parent = parent

    @property_cached
    def count(self):
        """Check both read counts are equal and return that number."""
        assert self.fwd.count == self.rev.count
        return self.fwd.count

    def open(self):
        self.fwd.open()
        self.rev.open()

    def parse(self):
        return izip(self.fwd.parse(), self.rev.parse())

    def close(self):
        self.fwd.close()
        self.rev.close()

    def create(self):
        self.fwd.create()
        self.rev.create()
        return self

    def add(self, f, r):
        return self.add_pair((f,r))

    def add_pair(self, pair):
        self.fwd.add_seq(pair[0])
        self.rev.add_seq(pair[1])

    def remove(self):
        self.fwd.remove()
        self.rev.remove()

    @property
    def progress(self):
        """Just like self.parse but display a progress bar."""
        return tqdm(self, total=len(self))

    def subsample(self, down_to, dest_pair=None):
        # Check size #
        assert down_to < len(self)
        # Make new pair of files #
        if dest_pair is None:
            dest_fwd_path = self.fwd.path.new_name_insert("subsampled")
            dest_rev_path = self.rev.path.new_name_insert("subsampled")
            dest_pair = self.__class__(dest_fwd_path, dest_rev_path)
        # Do it #
        dest_pair.create()
        for pair in isubsample(self, down_to): dest_pair.add_pair(pair)
        dest_pair.close()
        # Did it work #
        assert len(dest_pair) == down_to

    #------------------------------- Extensions ------------------------------#
    def parse_primers(self, *args, **kwargs):
        fwd_gen = self.fwd.parse_primers(*args, **kwargs)
        rev_gen = self.rev.parse_primers(*args, **kwargs)
        generator = izip(fwd_gen, rev_gen)
        return GenWithLength(generator, len(fwd_gen))

###############################################################################
class PairedFASTQ(PairedFASTA):
    """Read and write FASTQ file pairs without using too much RAM."""
    format = 'fastq'

    def __init__(self, fwd, rev, parent=None):
        # FASTQ objects #
        self.fwd = FASTQ(fwd)
        self.rev = FASTQ(rev)
        # Extra #
        self.gzipped = self.fwd.gzipped
        self.parent = parent

    def validate(self):
        """Call fastQValidator on these files."""
        self.fwd.validator()
        self.rev.validator()
