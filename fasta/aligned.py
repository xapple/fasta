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

    def gblocks(self,
                new_path = None,
                seq_type = 'nucl' or 'prot'):
        """Apply the gblocks filtering algorithm to the alignment.
        See http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html"""
        # Run it #
        if seq_type == 'nucl': t_option = "-t=d"
        if seq_type == 'prot': t_option = "-t=p"
        # Run it #
        sh.gblocks91(self, t_option, '-p=n', "-b5=a", _ok_code=[0,1])
        created_file = self.path + '-gb'
        assert os.path.exists(created_file)
        # Replace it maybe #
        if new_path is None: new_path = self.path
        # Reformat FASTA to the new path #
        AlignIO.write(AlignIO.parse(created_file, 'fasta'), new_path, 'fasta')
        # Clean up #
        os.remove(created_file)

    def build_tree(self,
                   new_path    = None,
                   seq_type    = 'nucl' or 'prot',
                   num_threads = None,
                   free_cores  = 2,
                   keep_dir    = False):
        """Make a tree with raxml. Note that you need at least four
        taxa to express some evolutionary history on an unrooted tree"""
        # Check length #
        assert len(self) > 3
        # Check output #
        if new_path is None: new_path = self.prefix_path + '.tree'
        # What model to choose #
        if seq_type == 'nucl': model = "GTRGAMMA"
        if seq_type == 'prot': model = "PROTGAMMAJTTF"
        # Threads #
        if num_threads is None: num_threads = multiprocessing.cpu_count() - free_cores
        else:                   num_threads = int(num_threads) - free_cores
        num_threads = max(1, num_threads)
        # Run it #
        temp_dir = new_temp_dir()
        sh.raxml811('-m', model, "-T", num_threads, '-p', 1, '-s', self.path, '-n', 'tree', '-w', temp_dir, '-f', 'a', '-x', 1, '-N', 'autoMR')
        # Move into place #
        if keep_dir: shutil.move(temp_dir, new_path)
        else:        shutil.move(temp_dir + 'RAxML_bestTree.tree', new_path)
        # Return #
        return FilePath(new_path)