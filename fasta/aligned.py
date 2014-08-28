# Built-in modules #
import os, multiprocessing
from collections import OrderedDict

# Internal modules #
from fasta import FASTA
from plumbing.cache import property_cached
from plumbing.autopaths import FilePath
from plumbing.tmpstuff import new_temp_dir

# Third party modules #
import sh, shutil
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

################################################################################
class AlignedFASTA(FASTA):
    """Also a FASTA file technically, but contains an alignment"""
    extension = 'aln'

    def __iter__(self): return iter(self.parse())

    def parse(self):
        self.open()
        return AlignIO.read(self.handle, self.format)

    def add_seq(self, seq):
        self.buffer.append(seq)

    def flush(self):
        AlignIO.write(MultipleSeqAlignment(self.buffer), self.handle, self.format)

    def write(self, reads):
        if not self.directory.exists: self.directory.create()
        self.open('w')
        AlignIO.write(reads, self.handle, self.format)
        self.close()

    @property_cached
    def sequences(self):
        return OrderedDict(((seq.id, seq) for seq in self))

    def gblocks(self, new_path=None, nucleotide=True):
        """Apply the gblocks filtering algorithm to the alignment and overwrite it.
        See http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html"""
        # Run it #
        if nucleotide: t = "-t=d"
        else: t = "-t=p"
        # Run it #
        sh.gblocks91(self, t, '-p=n', "-b5=a", _ok_code=[0,1])
        created_file = self.path + '-gb'
        assert os.path.exists(created_file)
        # Replace it #
        if new_path is None:
            os.remove(self.path)
            new_path = self.path
        # Reformat FASTA #
        AlignIO.write(AlignIO.parse(created_file, 'fasta'), new_path, 'fasta')
        # Clean up #
        os.remove(created_file)

    def build_tree(self, out_path=None, nucleotide=True):
        """Make a tree with raxml. Note that you need at least four
        taxa to express some evolutionary history on an unrooted tree"""
        # Check length #
        assert len(self) > 3
        # Check output #
        if out_path is None: out_path = self.prefix_path + '.tree'
        elif not isinstance(out_path, FilePath): out_path = FilePath(out_path)
        # Run it #
        temp_dir = new_temp_dir()
        cpus = max(multiprocessing.cpu_count(), 4) - 2
        if nucleotide: model = "GTRGAMMA"
        else: model = "PROTGAMMAJTTF"
        sh.raxml811('-m', model, "-T", cpus, '-p', 1, '-s', self.path, '-n', 'tree', '-w', temp_dir)
        # Move into place #
        shutil.move(temp_dir + 'RAxML_parsimonyTree.tree', out_path)