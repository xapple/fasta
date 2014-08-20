# Built-in modules #
import os, shutil

# Internal modules #
from fasta import FASTA
from plumbing.common import imean
from plumbing.cache import property_cached
from plumbing.autopaths import FilePath, DirectoryPath
from plumbing.tmpstuff import new_temp_dir

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
        if self.gziped: return int(sh.zgrep('-c', "^+$", self.path, _ok_code=[0,1]))
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
        mean = imean(s for r in self for s in r.letter_annotations["phred_quality"])
        self.close()
        return mean

    def fastqc(self, directory=None):
        # Default case #
        if directory is None:
            sh.fastqc(self.path, '-q')
            os.remove(self.prefix_path + '_fastqc.zip')
            return DirectoryPath(self.prefix.split('.')[0] + '_fastqc/')
        # Case directory #
        if directory is not None:
            if not isinstance(directory, DirectoryPath): directory = DirectoryPath(directory)
            if directory.exists: directory.remove()
            tmp_dir = new_temp_dir()
            sh.fastqc(self.path, '-q', '-o', tmp_dir)
            created_dir = tmp_dir + self.prefix.split('.')[0] + '_fastqc/'
            shutil.move(created_dir, directory)
            tmp_dir.remove()
            return directory

