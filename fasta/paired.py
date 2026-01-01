#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, itertools
from functools import cached_property

# Internal modules #
from fasta import FASTA, FASTQ

# First party modules #
from autopaths.file_size import FileSize
from autopaths.dir_path import DirectoryPath
from plumbing.common import isubsample, GenWithLength

# Third party modules #
from tqdm import tqdm

###############################################################################
class PairedFASTA:
    """Read and write FASTA file pairs without using too much RAM."""

    format = 'fasta'
    base_class = FASTA

    def __init__(self, fwd, rev, parent=None):
        # FASTA objects #
        self.fwd = self.base_class(fwd)
        self.rev = self.base_class(rev)
        # Set the prefixes #
        self.fwd.pair = 'fwd'
        self.rev.pair = 'rev'
        # Extra #
        self.parent = parent

    def __len__(self):  return self.count
    def __iter__(self): return self.parse()
    def __bool__(self): return bool(self.fwd) and bool(self.rev)
    def __repr__(self): return '<%s object on "%s" and "%s">' % \
                        (self.__class__.__name__, self.fwd.path, self.rev.path)

    def __enter__(self): return self.create()
    def __exit__(self, exc_type, exc_value, traceback): self.close()

    #------------------------------ Properties -------------------------------#
    @property
    def exists(self):
        return self.fwd.exists and self.rev.exists

    @property
    def prefix_path(self):
        return os.path.commonprefix((self.fwd.path, self.rev.path))

    @property
    def gzipped(self):
        assert self.fwd.gzipped == self.rev.gzipped
        return self.fwd.gzipped

    @cached_property
    def count(self):
        return self.fwd.count

    @property
    def progress(self):
        """Just like `self.parse` but display a progress bar."""
        return tqdm(self, total=len(self))

    @property
    def size(self):
        """Human-readable file size."""
        return FileSize(self.fwd.count_bytes + self.rev.count_bytes)

    @property
    def lengths(self):
        """All the sequence lengths, one by one, in an iterator."""
        return map(next, itertools.cycle((self.fwd.lengths, self.rev.lengths)))

    @property
    def first(self):
        """Just the first sequence pair."""
        return self.fwd.first, self.rev.first

    # Steal the `graphs` property from the base class (partial inheritance)
    graphs = FASTA.graphs

    @property
    def count_cache_path(self):
        return self.fwd.count_cache_path

    @count_cache_path.setter
    def count_cache_path(self, new_path):
        # Let's just count the forward reads instead of counting both files #
        self.fwd.count_cache_path = new_path
        # If we have a directory, then we can cache all the files #
        if isinstance(new_path, DirectoryPath):
            self.rev.count_cache_path = new_path
            self.singletons.count_cache_path = new_path
            self.other_reads.count_cache_path = new_path

    #------------------------------- Methods ---------------------------------#
    def open(self):
        self.fwd.open()
        self.rev.open()

    def parse(self):
        return zip(self.fwd.parse(), self.rev.parse())

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

    def compress(self, *args, **kwargs):
        # Make two separate calls #
        self.fwd.compress(*args, **kwargs)
        self.rev.compress(*args, **kwargs)
        # Print the total size #
        if kwargs.get('verbose'): print("Total size: %s." % self.size)

    def check_counts_match(self):
        """Check that both read counts are equal and return that number."""
        assert self.fwd.count == self.rev.count
        return self.fwd.count

    def get_id(self, fwd_id=None, rev_id=None, progress=False):
        """
        Extract one sequence from the file based on its ID.
        This is highly ineffective on large files.
        Consider using the SQLite API instead or memory map the file.
        """
        # Display a progress bar #
        if progress is True:  progress = tqdm
        if progress is False: progress = lambda x:x
        # Check at least one ID is specified #
        if fwd_id is None and rev_id is None:
            msg = "Please specify at least a forward or reverse ID."
            raise Exception(msg)
        # Three cases: 1. Only Forward #
        if fwd_id is not None and rev_id is None:
            for fwd, rev in progress(self):
                if fwd_id == fwd.id: return fwd, rev
        # Three cases: 2. Only Reverse #
        if fwd_id is None and rev_id is not None:
            for fwd, rev in progress(self):
                if rev_id == rev.id: return fwd, rev
        # Three cases: 3. Both specified #
        if fwd_id is not None and rev_id is not None:
            for fwd, rev in progress(self):
                if fwd_id == fwd.id and rev_id == rev.id: return fwd, rev

    def subsample(self, down_to, dest_pair=None):
        # Check size #
        assert down_to < len(self)
        # Make a new pair of files #
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
        # Return #
        return dest_pair

    #------------------------------- Extensions ------------------------------#
    def parse_primers(self, *args, **kwargs):
        fwd_gen = self.fwd.parse_primers(*args, **kwargs)
        rev_gen = self.rev.parse_primers(*args, **kwargs)
        generator = zip(fwd_gen, rev_gen)
        return GenWithLength(generator, len(fwd_gen))

    #------------------------------ Extra files ------------------------------#
    def file_with_suffix(self, suffix=''):
        """
        Sometimes, an algorithm in a pipeline can produce singletons or extra
        reads with are not paired anymore. In such cases, we can automatically
        create a file to hold them in the same directory.
        """
        # If there is no commonality don't add a period #
        if not self.prefix_path.endswith('/') and suffix != '': dot = '.'
        else: dot = ''
        # Add the requested suffix #
        path = self.prefix_path + dot + suffix + '.' + self.format
        # Return #
        return self.base_class(path)

    @cached_property
    def singletons(self): return self.file_with_suffix('singletons')

    @cached_property
    def other_reads(self): return self.file_with_suffix('others')

###############################################################################
class PairedFASTQ(PairedFASTA):
    """Read and write FASTQ file pairs without using too much RAM."""

    format = 'fastq'
    base_class = FASTQ

    def validate(self):
        """Call fastQValidator on these files."""
        self.fwd.validator()
        self.rev.validator()