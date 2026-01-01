#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, sys, io, gzip, shutil, itertools
from collections import Counter, OrderedDict
from functools import cached_property

# First party modules #
from plumbing.common          import isubsample
from plumbing.color           import Color
from plumbing.cache           import property_pickled
from autopaths.file_path      import FilePath
from autopaths.tmp_path       import new_temp_path

# Third party modules #
import sh
from tqdm import tqdm

################################################################################
class FASTA(FilePath):
    """
    A single FASTA file somewhere in the filesystem. You can read from it in
    several convenient ways. You can write to it in a automatically buffered
    way. There are several other things you can do with a FASTA file.
    """

    format      = 'fasta'
    ext         = 'fasta'
    buffer_size = 1000
    debug       = False

    def __len__(self): return self.count

    def __repr__(self):
        return '<%s object on "%s">' % (self.__class__.__name__, self.path)

    def __contains__(self, other): return other in self.ids

    def __enter__(self): return self.create()

    def __exit__(self, exc_type, exc_value, traceback): self.close()

    def __iter__(self):
        for seq in self.parse(): yield seq
        self.close()

    def __getitem__(self, key):
        if   isinstance(key, str): return self.sequences[key]
        elif isinstance(key, int): return self.sequences.items()[key]
        elif isinstance(key, slice):
            return itertools.islice(self, key.start, key.stop, key.step)

    #----------------------------- Properties --------------------------------#
    @property
    def gzipped(self): return True if self.path.endswith('gz') else False

    @property
    def first(self):
        """Just the first sequence."""
        from Bio import SeqIO
        self.open()
        seq = next(SeqIO.parse(self.handle, self.format))
        self.close()
        return seq

    count_cache_path = None
    @cached_property
    def count(self):
        """
        Count the total number sequences in the current FASTA files.
        For some files this can take some time, so you can set the attribute
        `self.count_cache_path` to point to a pickle file in which the result
        will be memoized.
        Note: we should probably check for file size changes instead of just
        caching once.
        """
        # In case we want to record the result to disk #
        if self.count_cache_path is not None:
            memoize = property_pickled(self.__class__.count_total,
                                       path=self.count_cache_path)
            return memoize.__get__(self, self.__class__)
        # Default case #
        return self.count_total()

    def count_total(self):
        # Message for debugging purposes #
        if self.debug:
            print("-> counting reads in `%s`" % self.path)
        # If we are gzipped we can just use zgrep #
        if self.gzipped:
            return int(sh.zgrep('-c', "^>", self.path, _ok_code=[0,1]))
        else:
            return int(sh.grep('-c', "^>", self.path, _ok_code=[0,1]))

    @property
    def lengths(self):
        """All the sequence lengths, one by one, in an iterator."""
        return map(len, self)

    @cached_property
    def lengths_counter(self):
        """A Counter() object with all the lengths inside."""
        return Counter((len(s) for s in self.parse()))

    #-------------------------- Basic IO methods -----------------------------#
    def open(self, mode='r'):
        # Two cases #
        if self.gzipped:
            self.handle = gzip.open(self.path, mode)
            self.handle = io.TextIOWrapper(self.handle, encoding='utf8')
        else:
            self.handle = open(self.path, mode)
        # For convenience #
        return self.handle

    def close(self):
        # Case we were writing to the file #
        if hasattr(self, 'buffer'):
            self.flush()
            del self.buffer
        # Standard case #
        self.handle.close()
        # For pickling purposes (can't use dill on gzip handles) #
        del self.handle

    def parse(self):
        self.open()
        from Bio import SeqIO
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
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        self.add_seq(SeqRecord(Seq(seq), id=name, description=description))

    def add_fasta(self, path):
        """Use this method to add an other fasta to this fasta."""
        path = FASTA(path)
        self.add(path)

    def add_fastas(self, paths):
        """Use this method to add a bunch of fastas to this fasta."""
        for p in paths: self.add_fasta(p)

    def flush(self):
        """Empty the buffer."""
        from Bio import SeqIO
        for seq in self.buffer:
            SeqIO.write(seq, self.handle, self.format)
        self.buffer = []

    def write(self, reads):
        from Bio import SeqIO
        if not self.directory.exists: self.directory.create()
        self.open('w')
        SeqIO.write(reads, self.handle, self.format)
        self.close()
        return self

    #-------------------------- Compressing the data -------------------------#
    def compress(self,
                 new_path    = None,
                 remove_orig = False,
                 method      = 'pigz',
                 verbose     = False):
        """Turn this FASTA file into a gzipped FASTA file."""
        # Check we are not compressed already #
        if self.gzipped and new_path is not False:
            msg = "The %s file '%s' is already compressed."
            raise Exception(msg % (self.format, self.path))
        # Print a message #
        if verbose:
            msg = "Compressing to '%s' with method '%s'."
            print(msg % (self.path, method))
        # Call method #
        return self.gzip_to(new_path, remove_orig, method)

    #------------------------- When IDs are important ------------------------#
    @cached_property
    def ids(self):
        """A frozen set of all unique IDs in the file."""
        as_list = [seq.description.split()[0] for seq in self]
        as_set = frozenset(as_list)
        assert len(as_set) == len(as_list)
        return as_set

    def get_id(self, id_num):
        """
        Extract one sequence from the file based on its ID.
        This is highly ineffective on large files.
        Consider using the SQLite API instead or memory map the file.
        """
        for seq in self:
            if seq.id == id_num: return seq

    @cached_property
    def sequences(self):
        """
        Another way of easily retrieving sequences. Also highly ineffective.
        Consider using the SQLite API instead.
        """
        return OrderedDict(((seq.id, seq) for seq in self))

    @cached_property
    def sql(self):
        """
        If you access this attribute, we will build an SQLite database
        out of the FASTA file. You will be able to access everything in an
        indexed fashion, and use the blaze library via `sql.frame`.
        """
        from fasta.indexed import DatabaseFASTA, fasta_to_sql
        db = DatabaseFASTA(self.prefix_path + ".db")
        if not db.exists: fasta_to_sql(self.path, db.path)
        return db

    @cached_property
    def length_by_id(self):
        """
        In some use cases you just need the sequence lengths in an indexed
        fashion. If you access this attribute, we will make a hash map in
        memory.
        """
        hash_map = dict((seq.id, len(seq)) for seq in self)
        tmp = hash_map.copy()
        hash_map.update(tmp)
        return hash_map

    #----------------- Ways of interacting with the data --------------------#
    def subsample(self, down_to=1, new_path=None, verbose=True):
        """Pick a given number of sequences from the file pseudo-randomly."""
        # Pick the destination path #
        if new_path is None:
            subsampled = self.__class__(new_temp_path())
        elif isinstance(new_path, FASTA):
            subsampled = new_path
        else:
            subsampled = self.__class__(new_path)
        # Check size #
        if down_to > len(self):
            message = "Can't subsample %s down to %i. Only down to %i."
            print(Color.ylw + message % (self, down_to, len(self)) + Color.end)
            self.copy(new_path)
            return
        # Select verbosity #
        import tqdm
        if verbose: wrapper = lambda x: tqdm.tqdm(x, total=self.count)
        else: wrapper = lambda x: x
        # Generator #
        def iterator():
            for read in wrapper(isubsample(self, down_to)):
                yield read
        # Do it #
        subsampled.write(iterator())
        # Did it work #
        assert len(subsampled) == down_to
        # Return #
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
        # Replace it #
        if new_path is None:
            os.remove(self.path)
            shutil.move(numbered, self.path)
        # Return #
        return numbered

    def rename_with_prefix(self, prefix="", new_path=None, in_place=True,
                           remove_desc=True):
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
        # Replace it #
        if in_place:
            os.remove(self.path)
            shutil.move(prefixed, self.path)
        # Return #
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

    def extract_length(self, lower_bound=None,
                             upper_bound=None,
                             new_path=None):
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
        # Return #
        return fraction

    def extract_sequences(self, ids,
                          new_path = None,
                          in_place = False,
                          verbose  = False):
        """
        Will take all the sequences from the current file who's id appears in
        the ids given and place them in a new file.
        If no path is given, a new temporary path is created and returned.
        If `in_place` is set to True, the original file is removed and replaced
        with the result of the extraction.
        Optionally, the argument `ids` can be a function which has to take
        one string as only input and return True for keeping the sequence and
        False for discarding the sequence.
        """
        # Temporary path #
        if new_path is None: new_fasta = self.__class__(new_temp_path())
        elif isinstance(new_path, FASTA): new_fasta = new_path
        else:                new_fasta = self.__class__(new_path)
        # Select verbosity #
        import tqdm
        wrapper = tqdm.tqdm if verbose else lambda x: x
        # Simple generator #
        def simple_match(reads):
            for r in wrapper(reads):
                if r.id in ids: yield r
        # Generator with function #
        def function_match(reads):
            for r in wrapper(reads):
                if ids(r.id): yield r
        # Do it #
        if callable(ids):
            new_fasta.write(function_match(self))
        else:
            new_fasta.write(simple_match(self))
        # Return #
        if in_place:
            os.remove(self.path)
            shutil.move(new_fasta, self.path)
            return self
        else: return new_fasta

    def remove_trailing_stars(self, new_path=None, in_place=True, check=False):
        """
        Remove the bad character that can be inserted by some programs at the
        end of sequences.
        """
        # Optional check #
        if check and int(sh.grep('-c', '\\*', self.path, _ok_code=[0,1])) == 0:
            return self
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
        # Return #
        return new_fasta

    def _generator_mod(self, generator, new_path=None, in_place=True):
        """
        Generic way of modifying the current fasta either in place or
        with a new destination pass.
        Simply, pass a generator function that will yield the new sequences
        given the current ones.
        """
        # Temporary path #
        if new_path is None: new_fasta = self.__class__(new_temp_path())
        elif isinstance(new_path, FASTA): new_fasta = new_path
        else: new_fasta = self.__class__(new_path)
        # Do it #
        new_fasta.write(generator())
        # Return #
        if in_place:
            os.remove(self.path)
            shutil.move(new_fasta, self.path)
            return self
        else: return new_fasta

    def remove_duplicates(self, new_path=None, in_place=True):
        """
        If several entries have the same ID in the FASTA file, keep only the
        first appearance and remove all the others.
        """
        # Generator #
        def unique_entries():
            seen = set()
            for i, read in enumerate(self):
                if read.id in seen: continue
                else:
                    seen.add(read.id)
                    yield read
        # Return #
        return self._generator_mod(unique_entries, new_path, in_place)

    def convert_U_to_T(self, new_path=None, in_place=True):
        # Generator #
        def all_U_to_T():
            for i, read in enumerate(self):
                read.seq = read.seq.back_transcribe()
                yield read
        # Return #
        return self._generator_mod(all_U_to_T, new_path, in_place)

    #---------------------------- Third party programs -----------------------#
    def muscle_align(self, out_path=None):
        """We align the sequences in the fasta file with muscle."""
        if out_path is None: out_path = self.prefix_path + '.aln'
        sh.muscle38("-in", self.path, "-out", out_path)
        from fasta.aligned import AlignedFASTA
        return AlignedFASTA(out_path)

    def mothur_align(self, ref_path):
        """
        We align the sequences in the fasta file with mothur and a template.
        """
        # Run it #
        msg = "#align.seqs(candidate=%s, template=%s, search=blast," \
              "flip=false, processors=8);"
        sh.mothur(msg % (self.path, ref_path))
        # Move things #
        shutil.move(self.path[:-6] + '.align',        self.p.aligned)
        shutil.move(self.path[:-6] + '.align.report', self.p.report)
        shutil.move(self.path[:-6] + '.flip.accnos',  self.p.accnos)
        # Clean up #
        if os.path.exists('formatdb.log'):
            os.remove('formatdb.log')
        if os.path.exists('error.log') and os.path.getsize('error.log') == 0:
            os.remove('error.log')
        for path in sh.glob('mothur.*.logfile'):
            os.remove(path)
        # Return #
        return self.p.aligned

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

    def index_bwa(self, out_path=None, verbose=True):
        """
        Create an index on the fasta file compatible with BWA.
        This will take up approximately five times the disk space of the
        original gzipped FASTA file it is based on.
        """
        # By default, in the same directory #
        if out_path is None: out_path = self.path
        # Message #
        if verbose:
            msg = "Creating BWA index on '%s' at '%s'."
            print(msg % (self.path, out_path))
        # Command line options #
        options = {'p':    out_path,    # Prefix
                   'a':    'bwtsw'}     # Algorithm
        # Redirect output #
        if verbose:
            options['_out'] = sys.stdout
            options['_err'] = sys.stderr
        # Call the BWA executable #
        cmd = sh.bwa('index', self.path, **options)
        # Show the full command #
        if verbose: print("Ran the following command:\n  $ %s" % cmd.ran)
        # Return #
        return FilePath(out_path + '.bwt')

    #--------------------------------- Graphs --------------------------------#
    @cached_property
    def graphs(self):
        """
        Sorry for the black magic. The result is an object whose attributes
        are all the graphs found in `./graphs.py` initialized with this
        instance as only argument.
        """
        # Make a dummy object that is empty at first #
        class Graphs: pass
        result = Graphs()
        # Loop over graphs classes and add them as attributes #
        from fasta import graphs
        for graph in graphs.__all__:
            GraphClass = getattr(graphs, graph)
            setattr(result, GraphClass.short_name, GraphClass(self))
        # Return #
        return result

    #-------------------------------- Primers -------------------------------#
    def parse_primers(self, primers, mismatches=None):
        """
        Takes care of identifying primers inside every sequence.
        Instead of yielding Seq objects now we yield ReadWithPrimers objects.
        These have extra properties that show the start and end positions
        of all primers found.
        """
        # Default is zero #
        if mismatches is None: mismatches = 0
        # Get the search expressions with mismatches #
        from fasta.primers import PrimersRegexes
        regexes = PrimersRegexes(primers, mismatches)
        # Generate a new special object for every read #
        from fasta.primers import ReadWithPrimers
        read_with_primer = lambda read: ReadWithPrimers(read, regexes)
        generator = (read_with_primer(r) for r in self.parse())
        # Add the length to the generator #
        from plumbing.common import GenWithLength
        generator = GenWithLength(generator, len(self))
        # Return #
        return generator