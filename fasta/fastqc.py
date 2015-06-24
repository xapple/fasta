# Futures #
from __future__ import division

# Built-in modules #
import os, shutil

# First party modules #
from fasta import FASTQ
from plumbing.autopaths import DirectoryPath
from plumbing.tmpstuff import new_temp_dir
from plumbing.cache import property_cached
from plumbing.slurm import num_processors

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
        if dest is None: self.dest = DirectoryPath(self.source.prefix_path + '.fastqc')

    def run(self, cpus=None):
        # Check version #
        assert "v0.10.1" in sh.fastqc101('--version')
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Destination absent #
        if self.dest is None:
            sh.fastqc101(self.source, '-q', '-t', cpus)
            os.remove(self.source.prefix_path + '_fastqc.zip')
        # Destination given #
        if self.dest is not None:
            self.tmp_dir = new_temp_dir()
            sh.fastqc101(self.source, '-q', '-o', self.tmp_dir, '-t', cpus)
            created_name = self.source.prefix.split('.')[0] + '_fastqc/'
            created_dir  = self.tmp_dir + created_name
            if self.dest.exists: shutil.rmtree(self.dest)
            shutil.move(created_dir, self.dest)
            self.tmp_dir.remove()
        # Strange case where created_dir is not renamed to self.dest but put inside instead #
        if self.dest + created_name in self.dest.flat_directories:
            raise Exception("Looks like the shutil.move didn't do what was expected. Should fix.")
        # Return #
        return self.results

    @property
    def output_dir(self):
        if self.dest is None: return self.source.split('.')[0] + '_fastqc/'
        else: return self.dest

    @property_cached
    def results(self):
        results = FastQCResults(self.output_dir)
        if not results: raise Exception("You can't access results from FastQC before running the tool.")
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