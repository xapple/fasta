# Built-in modules #
from itertools import izip

# Internal modules #
from fasta import FASTA, FASTQ
from plumbing.common import isubsample
from plumbing.cache import property_cached

# Third party modules #
from tqdm import tqdm

###############################################################################
class PairedFASTA(object):
    """Read and write FASTA file pairs without using too much RAM"""

    def __len__(self): return self.count
    def __iter__(self): return self.parse()
    def __enter__(self): return self.create()
    def __exit__(self, exc_type, exc_value, traceback): self.close()
    def __repr__(self): return '<%s object on "%s" and "%s">' % \
                        (self.__class__.__name__, self.fwd.path, self.rev.path)

    @property
    def exists(self): return self.fwd.exists

    def __init__(self, fwd, rev, parent=None):
        # FASTA objects #
        self.fwd = FASTA(fwd)
        self.rev = FASTA(rev)
        # Extra #
        self.gzipped = self.fwd.gzipped
        self.parent = parent

    @property_cached
    def count(self):
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

    def add_pair(self, pair):
        self.fwd.add_seq(pair[0])
        self.rev.add_seq(pair[1])

    @property
    def progress(self):
        """Just like self.parse but display a progress bar"""
        return tqdm(self, total=len(self))

    def subsample(self, down_to, dest_pair=None):
        # Check size #
        assert down_to < len(self)
        # Make new pair of files #
        if dest_pair is None:
            dest_fwd_path = self.fwd_path.new_name_insert("subsampled")
            dest_rev_path = self.rev_path.new_name_insert("subsampled")
            dest_pair = self.__class__(dest_fwd_path, dest_rev_path)
        # Do it #
        dest_pair.create()
        for pair in isubsample(self, down_to): dest_pair.add_pair(pair)
        self.subsampled.close()
        # Did it work #
        assert len(dest_pair) == down_to

###############################################################################
class PairedFASTQ(PairedFASTA):
    """Read and write FASTQ file pairs without using too much RAM"""

    def __init__(self, fwd, rev, parent=None):
        # FASTQ objects #
        self.fwd = FASTQ(fwd)
        self.rev = FASTQ(rev)
        # Extra #
        self.gzipped = self.fwd.gzipped
        self.parent = parent
