#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, sys, gzip, shutil, itertools
from collections import Counter, OrderedDict
from six import string_types

# Internal modules #
from fasta import graphs, primers

# First party modules #
from plumbing.common     import isubsample
from plumbing.color      import Color
from plumbing.cache      import property_cached
from autopaths.file_path import FilePath
from autopaths.tmp_path  import new_temp_path

# Third party modules #
import sh
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

################################################################################
class FASTA(FilePath):
    """
    A single FASTA file somewhere in the filesystem. You can read from it in
    several convenient ways. You can write to it in a automatically buffered way.
    There are several other things you can do with a FASTA file. Look at the class.
    """

    format      = 'fasta'
    ext         = 'fasta'
    buffer_size = 1000

    def __len__(self): return self.count
    def __repr__(self): return '<%s object on "%s">' % (self.__class__.__name__, self.path)
    def __contains__(self, other): return other in self.ids

    def __enter__(self): return self.create()
    def __exit__(self, exc_type, exc_value, traceback): self.close()

    def __iter__(self):
        for seq in self.parse(): yield seq
        self.close()

    def __getitem__(self, key):
        if   isinstance(key, string_types): return self.sequences[key]
        elif isinstance(key, int):          return self.sequences.items()[key]
        elif isinstance(key, slice):        return itertools.islice(self, key.start, key.stop, key.step)

    #----------------------------- Properties --------------------------------#
    @property
    def gzipped(self): return True if self.path.endswith('gz') else False

    @property
    def first(self):
        """Just the first sequence"""
        self.open()
        seq = SeqIO.parse(self.handle, self.format).next()
        self.close()
        return seq

    @property_cached
    def count(self):
        """Should probably check for file size changes instead of just caching once #TODO"""
        if self.gzipped: return int(sh.zgrep('-c', "^>", self.path, _ok_code=[0,1]))
        else: return int(sh.grep('-c', "^>", self.path, _ok_code=[0,1]))

    @property
    def lengths(self):
        """All the lengths, one by one, in a list."""
        return map(len, self.parse())

    @property_cached
    def lengths_counter(self):
        """A Counter() object with all the lengths inside."""
        return Counter((len(s) for s in self.parse()))

    #-------------------------- Basic IO methods -----------------------------#
    def open(self, mode='r'):
        # Two cases #
        if self.gzipped: self.handle = gzip.open(self.path, mode)
        else:            self.handle =      open(self.path, mode)
        # For convenience #
        return self.handle

    def close(self):
        if hasattr(self, 'buffer'):
            self.flush()
            del self.buffer
        self.handle.close()

    def parse(self):
        self.open()
        return SeqIO.parse(self.handle, self.format)

    @property
    def progress(self):
        """Just like self.parse() but will display a progress bar."""
        return tqdm(self, total=len(self))

    def create(self):
        """Create the file on the file system."""
        self.buffer = []
        self.buf_count = 0
        if not self.directory.exists: self.directory.create()
        self.open('w')
        return self

    def add(self, seqs):
        """Use this method to add a bunch of SeqRecords at once."""
        for seq in seqs: self.add_seq(seq)

    def add_seq(self, seq):
        """Use this method to add a SeqRecord object to this fasta."""
        self.buffer.append(seq)
        self.buf_count += 1
        if self.buf_count % self.buffer_size == 0: self.flush()

    def add_str(self, seq, name=None, description=""):
        """Use this method to add a sequence as a string to this fasta."""
        self.add_seq(SeqRecord(Seq(seq), id=name, description=description))

    def add_fasta(self, path):
        """Use this method to add an other fasta to this fasta."""
        path = FASTA(path)
        self.add(path)

    def add_fastas(self, paths):
        """Use this method to add bunch of fastas to this fasta."""
        for p in paths: self.add_fasta(p)

    def flush(self):
        """Empty the buffer."""
        for seq in self.buffer:
            SeqIO.write(seq, self.handle, self.format)
        self.buffer = []

    def write(self, reads):
        if not self.directory.exists: self.directory.create()
        self.open('w')
        SeqIO.write(reads, self.handle, self.format)
        self.close()

    #------------------------- When IDs are important ------------------------#
    @property_cached
    def ids(self):
        """A frozen set of all unique IDs in the file."""
        as_list = [seq.description.split()[0] for seq in self]
        as_set = frozenset(as_list)
        assert len(as_set) == len(as_list)
        return as_set

    def get_id(self, id_num):
        """
        Extract one sequence from the file based on its ID. This is highly ineffective.
        Consider using the SQLite API instead or memory map the file.
        """
        for seq in self:
            if seq.id == id_num: return seq

    @property_cached
    def sequences(self):
        """
        Another way of easily retrieving sequences. Also highly ineffective.
        Consider using the SQLite API instead.
        """
        return OrderedDict(((seq.id, seq) for seq in self))

    @property_cached
    def sql(self):
        """
        If you access this attribute, we will build an SQLite database
        out of the FASTA file and you will be able access everything in an
        indexed fashion, and use the blaze library via sql.frame
        """
        from fasta.indexed import DatabaseFASTA, fasta_to_sql
        db = DatabaseFASTA(self.prefix_path + ".db")
        if not db.exists: fasta_to_sql(self.path, db.path)
        return db

    @property_cached
    def length_by_id(self):
        """
        In some use cases you just need the sequence lengths in an indexed
        fashion. If you access this attribute, we will make a hash map in memory.
        """
        hash_map = dict((seq.id, len(seq)) for seq in self)
        tmp = hash_map.copy()
        hash_map.update(tmp)
        return hash_map

    #----------------- Ways of interacting with the data --------------------#
    def subsample(self, down_to=1, new_path=None):
        """Pick a number of sequences from the file pseudo-randomly."""
        # Auto path #
        if new_path is None: subsampled = self.__class__(new_temp_path())
        elif isinstance(new_path, FASTA): subsampled = new_path
        else:                subsampled = self.__class__(new_path)
        # Check size #
        if down_to > len(self):
            message = "Can't subsample %s down to %i. Only down to %i."
            print(Color.ylw + message % (self, down_to, len(self)) + Color.end)
            self.copy(new_path)
            return
        # Do it #
        subsampled.create()
        for seq in isubsample(self, down_to): subsampled.add_seq(seq)
        subsampled.close()
        # Did it work #
        assert len(subsampled) == down_to
        return subsampled

    def rename_with_num(self, prefix="", new_path=None, remove_desc=True):
        """Rename every sequence based on a prefix and a number."""
        # Temporary path #
        if new_path is None: numbered = self.__class__(new_temp_path())
        else:                numbered = self.__class__(new_path)
        # Generator #
        def numbered_iterator():
            for i,read in enumerate(self):
                read.id  = prefix + str(i)
                read.seq = read.seq.upper()
                if remove_desc: read.description = ""
                yield read
        # Do it #
        numbered.write(numbered_iterator())
        numbered.close()
        # Replace it #
        if new_path is None:
            os.remove(self.path)
            shutil.move(numbered, self.path)
        return numbered

    def rename_with_prefix(self, prefix="", new_path=None, in_place=True, remove_desc=True):
        """Rename every sequence based on a prefix."""
        # Temporary path #
        if new_path is None: prefixed = self.__class__(new_temp_path())
        else:                prefixed = self.__class__(new_path)
        # Generator #
        def prefixed_iterator():
            for i,read in enumerate(self):
                read.id = prefix + read.id
                if remove_desc: read.description = ""
                yield read
        # Do it #
        prefixed.write(prefixed_iterator())
        prefixed.close()
        # Replace it #
        if in_place:
            os.remove(self.path)
            shutil.move(prefixed, self.path)
        return prefixed

    def rename_sequences(self, mapping, new_path=None, in_place=False):
        """
        Will rename all sequences in the current fasta file using
        the mapping dictionary also provided. In place or at a new path.
        """
        # Where is the new file #
        if new_path is None: new_fasta = self.__class__(new_temp_path())
        else:                new_fasta = self.__class__(new_path)
        # Do it #
        new_fasta.create()
        for seq in self:
            new_name = mapping[seq.description]
            nucleotides = str(seq.seq)
            new_fasta.add_str(nucleotides, new_name)
        new_fasta.close()
        # Return #
        if in_place:
            os.remove(self.path)
            shutil.move(new_fasta, self.path)
            return self
        else: return new_fasta

    def extract_length(self, lower_bound=None, upper_bound=None, new_path=None):
        """Extract a certain length fraction and place them in a new file."""
        # Temporary path #
        if new_path is None: fraction = self.__class__(new_temp_path())
        elif isinstance(new_path, FASTA): fraction = new_path
        else:                fraction = self.__class__(new_path)
        # Generator #
        if lower_bound is None: lower_bound = 0
        if upper_bound is None: upper_bound = sys.maxsize
        def fraction_iterator():
            for read in self:
                if lower_bound <= len(read) <= upper_bound:
                    yield read
        # Do it #
        fraction.write(fraction_iterator())
        fraction.close()
        return fraction

    def extract_sequences(self, ids, new_path=None):
        """
        Will take all the sequences from the current file who's id appears in
        the ids given and place them in the new file path given.
        """
        # Temporary path #
        if new_path is None: new_fasta = self.__class__(new_temp_path())
        elif isinstance(new_path, FASTA): new_fasta = new_path
        else:                new_fasta = self.__class__(new_path)
        # Do it #
        new_fasta.create()
        for seq in self:
            if seq.id in ids: new_fasta.add_seq(seq)
        new_fasta.close()
        return new_fasta

    def remove_trailing_stars(self, new_path=None, in_place=True, check=False):
        """
        Remove the bad character that can be inserted by some programs at the
        end of sequences.
        """
        # Optional check #
        if check and int(sh.grep('-c', '\\*', self.path, _ok_code=[0,1])) == 0: return self
        # Faster with bash utilities #
        if in_place is True:
            sh.sed('-i', 's/\\*$//g', self.path)
            return self
        # Standard way #
        if new_path is None: new_fasta = self.__class__(new_temp_path())
        else:                new_fasta = self.__class__(new_path)
        new_fasta.create()
        for seq in self: new_fasta.add_str(str(seq.seq).rstrip('*'), seq.id)
        new_fasta.close()
        return new_fasta

    #---------------------------- Third party programs -----------------------#
    def align(self, out_path=None):
        """We align the sequences in the fasta file with muscle."""
        if out_path is None: out_path = self.prefix_path + '.aln'
        sh.muscle38("-in", self.path, "-out", out_path)
        return AlignedFASTA(out_path)

    def template_align(self, ref_path):
        """We align the sequences in the fasta file with mothur and a template."""
        # Run it #
        sh.mothur("#align.seqs(candidate=%s, template=%s, search=blast, flip=false, processors=8);" % (self.path, ref_path))
        # Move things #
        shutil.move(self.path[:-6] + '.align',        self.p.aligned)
        shutil.move(self.path[:-6] + '.align.report', self.p.report)
        shutil.move(self.path[:-6] + '.flip.accnos',  self.p.accnos)
        # Clean up #
        if os.path.exists('formatdb.log'): os.remove('formatdb.log')
        if os.path.exists('error.log') and os.path.getsize('error.log') == 0: os.remove('error.log')
        for p in sh.glob('mothur.*.logfile'): os.remove(p)

    def index_bowtie(self):
        """Create an index on the fasta file compatible with bowtie2."""
        # It returns exit code 1 if the fasta is empty #
        assert self
        # Call the bowtie executable #
        sh.bowtie2_build(self.path, self.path)
        return FilePath(self.path + '.1.bt2')

    def index_samtools(self):
        """Create an index on the fasta file compatible with samtools."""
        sh.samtools('faidx', self.path)
        return FilePath(self.path + '.fai')

    #--------------------------------- Graphs --------------------------------#
    @property_cached
    def graphs(self):
        """
        Sorry for the black magic. The result is an object whose attributes
        are all the graphs found in graphs.py initialized with this instance as
        only argument.
        """
        class Dummy: pass
        result = Dummy()
        for graph in graphs.__all__:
            cls = getattr(graphs, graph)
            setattr(result, cls.short_name, cls(self))
        return result

    @property_cached
    def length_dist(self):
        graph = self.graphs.length_dist
        if not graph: graph.plot()
        return graph

    #------------------------------- Extensions ------------------------------#
    def parse_primers(self, *args, **kwargs):
        return primers.parse_primers(self, *args, **kwargs)