# Built-in modules #

# Internal modules #
from fasta import FASTA
from plumbing.common import average
from plumbing.cache import property_cached
from plumbing.autopaths import FilePath

# Third party modules #
import sh
from Bio import SeqIO

################################################################################
class FASTQ(FASTA):
    """A single FASTQ file somewhere in the filesystem"""

    extension = 'fastq'
    format    = 'fastq'

    @property_cached
    def count(self):
        if self.gzipped: return int(sh.zgrep('-c', "^+$", self.path, _ok_code=[0,1]))
        return int(sh.grep('-c', "^+$", self.path, _ok_code=[0,1]))

    def to_fasta(self, path):
        with open(path, 'w') as handle:
            for r in self: SeqIO.write(r, handle, 'fasta')
        return FASTA(path)

    def to_qual(self, path):
        with open(path, 'w') as handle:
            for r in self: SeqIO.write(r, handle, 'qual')
        return FilePath(path)

    @property_cached
    def avg_quality(self):
        mean = average(s for r in self for s in r.letter_annotations["phred_quality"])
        self.close()
        return mean

    @property_cached
    def fastqc(self):
        from fasta.fastqc import FastQC
        return FastQC(self)