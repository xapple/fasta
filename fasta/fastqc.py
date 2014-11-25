# Futures #
from __future__ import division

# Built-in modules #
import os, shutil

# Internal modules #
from fasta import FASTQ
from plumbing.autopaths import DirectoryPath
from plumbing.tmpstuff import new_temp_dir
from plumbing.cache import property_cached

# Third party modules #
import sh

###############################################################################
class FastQC(object):
    """Takes care of running the FastQC program.
    See http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    Expects version 0.10.1."""

    def __repr__(self): return '<%s object on "%s">' % (self.__class__.__name__, self.source.path)

    def __init__(self, source, dest=None):
        # Basic #
        self.source = FASTQ(source)
        self.dest = DirectoryPath(dest)
        # Default case #
        if dest is None:
            self.dest = DirectoryPath(self.source.prefix_path + '.fastqc')

    def check(self):
        assert sh.fastqc('--v', )

    def run(self):
        if self.dest is None:
            sh.fastqc101(self.source, '-q')
            os.remove(self.source.prefix_path + '_fastqc.zip')
        if self.dest is not None:
            if self.dest.exists: self.dest.remove()
            self.tmp_dir = new_temp_dir()
            sh.fastqc101(self.source, '-q', '-o', self.tmp_dir)
            created_dir = self.tmp_dir + self.source.prefix.split('.')[0] + '_fastqc/'
            shutil.move(created_dir, self.dest)
            self.tmp_dir.remove()
            return self.results

    @property
    def output_dir(self):
        if self.dest is None: return self.source.split('.')[0] + '_fastqc/'
        else: return self.dest

    @property_cached
    def results(self):
        results = FastQCResults(self.output_dir)
        if not results: self.run()
        return results

###############################################################################
class FastQCResults(DirectoryPath):
    """A directory with the results from FastQC"""

    all_paths = """
    /Images/per_base_quality.png
    /Images/per_sequence_quality.png
    """

    def __nonzero__(self): return self.per_base_qual.exists

    @property
    def per_base_qual(self): return self.p.per_base_quality
    @property
    def per_seq_qual(self): return self.p.per_sequence_quality