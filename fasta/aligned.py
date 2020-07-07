#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, multiprocessing
from collections import OrderedDict

# Internal modules #
from fasta import FASTA
from plumbing.cache import property_cached
from autopaths.file_path import FilePath
from autopaths.tmp_path  import new_temp_path, new_temp_dir

# Third party modules #
import sh, shutil
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

################################################################################
class AlignedFASTA(FASTA):
    """Also a FASTA file technically, but contains an alignment."""
    ext = 'aln'

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
        See http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html
        Need to rename all sequences, because it will complain with long names."""
        # Temporary path #
        if new_path is None: final = self.__class__(new_temp_path())
        else:                final = self.__class__(new_path)
        # Mapping every sequence name with a random name #
        orig_name_to_temp = {seq.description: 'name' + str(i) for i,seq in enumerate(self)}
        temp_name_to_orig = {v: k for k, v in orig_name_to_temp.items()}
        # Rename every sequence with a random name #
        temp_fasta = self.rename_sequences(orig_name_to_temp)
        # Options #
        if seq_type == 'nucl': t_option = "-t=d"
        if seq_type == 'prot': t_option = "-t=p"
        # Run it #
        result = sh.gblocks91(temp_fasta.path, t_option, '-p=n', "-b4=3", "-b3=20", "-b5=a", _ok_code=[0,1])
        created_file = temp_fasta.path + '-gb'
        assert os.path.exists(created_file)
        # Check errors #
        if "Execution terminated" in result.stdout: raise Exception("gblocks crashed again.")
        # Back #
        temp_fasta.rename_sequences(temp_name_to_orig, final)
        # Return #
        return final

    def build_tree(self, *args, **kwargs):
        """Dispatch a tree build call. Note that you need at least four
        taxa to express some evolutionary history on an unrooted tree."""
        # Check length #
        assert len(self) > 3
        # Default option #
        algorithm = kwargs.pop(kwargs, None)
        if algorithm is None: algorithm = 'raxml'
        # Dispatch #
        if algorithm == 'raxml':    return self.build_tree_raxml(*args, **kwargs)
        if algorithm == 'fasttree': return self.build_tree_fast(*args, **kwargs)

    def build_tree_raxml(self,
                   new_path    = None,
                   seq_type    = 'nucl' or 'prot',
                   num_threads = None,
                   free_cores  = 2,
                   keep_dir    = False):
        """Make a tree with RAxML."""
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
        if keep_dir:
            shutil.rmtree(new_path)
            shutil.move(temp_dir, new_path)
        if not keep_dir:
            shutil.move(temp_dir + 'RAxML_bestTree.tree', new_path)
        # Return #
        return FilePath(new_path)

    def build_tree_fast(self, new_path=None, seq_type='nucl' or 'prot'):
        """Make a tree with FastTree. Names will be truncated however."""
        # Check output #
        if new_path is None: new_path = self.prefix_path + '.tree'
        # Command #
        command_args = []
        if seq_type == 'nucl': command_args += ['-nt']
        command_args += ['-gamma']
        command_args += ['-out', new_path]
        command_args += [self.path]
        # Run it #
        sh.FastTree(*command_args)
        # Return #
        return FilePath(new_path)

