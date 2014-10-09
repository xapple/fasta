b'This module needs Python 2.7.x'

# Special variables #
__version__ = '1.0.2'

# Built-in modules #
import os, gzip, shutil
from collections import Counter, OrderedDict

# Internal modules #
from fasta import graphs
from plumbing.common import isubsample
from plumbing.color import Color
from plumbing.cache import property_cached
from plumbing.autopaths import FilePath
from plumbing.tmpstuff import new_temp_path

# Third party modules #
import sh
from tqdm import tqdm
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

################################################################################
class FASTA(FilePath):
    """A single FASTA file somewhere in the filesystem. You can read from it in
    several convenient ways. You can write to it in a automatically buffered way.
    There are several other things you can do with a FASTA file. Look at the class."""

    format    = 'fasta'
    extension = 'fasta'
    buffer_size = 1000

    def __len__(self): return self.count
    def __iter__(self): return self.parse()
    def __repr__(self): return '<%s object on "%s">' % (self.__class__.__name__, self.path)
    def __nonzero__(self): return self.exists and os.path.getsize(self.path) != 0
    def __contains__(self, other): return other in self.ids

    def __enter__(self): return self.create()
    def __exit__(self, exc_type, exc_value, traceback): self.close()

    def __getitem__(self, key):
        if isinstance(key, basestring): return self.sequences[key]
        elif isinstance(key, int): return self.sequences.items()[key]

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
        """Should probably check file size instead of just caching once #TODO"""
        if self.gzipped: return int(sh.zgrep('-c', "^>", self.path, _ok_code=[0,1]))
        else: return int(sh.grep('-c', "^>", self.path, _ok_code=[0,1]))

    @property_cached
    def ids(self):
        return frozenset([seq.description.split()[0] for seq in self])

    @property
    def lengths(self):
        return map(len, self.parse())

    @property_cached
    def lengths_counter(self):
        return Counter((len(s) for s in self.parse()))

    def open(self, mode='r'):
        if self.gzipped: self.handle = gzip.open(self.path, mode)
        else: self.handle = open(self.path, mode)

    def close(self):
        if hasattr(self, 'buffer'):
            self.flush()
            del self.buffer
        self.handle.close()

    def parse(self):
        self.open()
        return SeqIO.parse(self.handle, self.format)

    def create(self):
        self.buffer = []
        self.buf_count = 0
        if not self.directory.exists: self.directory.create()
        self.open('w')
        return self

    def add_seq(self, seq):
        self.buffer.append(seq)
        self.buf_count += 1
        if self.buf_count % self.buffer_size == 0: self.flush()

    def add_str(self, seq, name=None, description=""):
        self.add_seq(SeqRecord(Seq(seq), id=name, description=description))

    def flush(self):
        for seq in self.buffer:
            SeqIO.write(seq, self.handle, self.format)
        self.buffer = []

    def write(self, reads):
        if not self.directory.exists: self.directory.create()
        self.open('w')
        SeqIO.write(reads, self.handle, self.format)
        self.close()

    @property
    def progress(self):
        """Just like self.parse but display a progress bar"""
        return tqdm(self, total=len(self))

    def get_id(self, id_num):
        """Extract one sequence from the file based on its ID. This is highly ineffective.
        Consider using the SQLite API instead."""
        for seq in self:
            if seq.id == id_num: return seq

    @property_cached
    def sequences(self):
        """Another way of easily retrieving sequences. Also highly ineffective.
        Consider using the SQLite API instead."""
        return OrderedDict(((seq.id, seq) for seq in self))

    @property_cached
    def sql(self):
        """If you access this attribute, we will build an SQLite database
        out of the FASTA file and you will be able access everything in an
        indexed fashion"""
        from fasta.database import DatabaseFASTA, fasta_to_sql
        db = DatabaseFASTA(self.prefix_path + ".db")
        if not db.exists: fasta_to_sql(self.path, db.path)
        return db

    @property_cached
    def length_by_id(self):
        """In some use cases you just need the sequence lengths in an indexed
        fashion. If you access this attribute, we will make a hash map in memory."""
        hashmap = dict((seq.id, len(seq)) for seq in self)
        tmp = hashmap.copy()
        hashmap.update(tmp)
        return hashmap

    def subsample(self, down_to=1, new_path=None):
        """Pick a number of sequences from the file randomly"""
        # Auto path #
        if not new_path: new_path = self.p.subsample
        # Check size #
        if down_to > len(self):
            message = "Can't subsample %s down to %i. Only down to %i."
            print Color.ylw + message % (self, down_to, len(self)) + Color.end
            self.copy(new_path)
            return
        # Make new file #
        self.subsampled = self.__class__(new_path)
        self.subsampled.create()
        # Do it #
        for seq in isubsample(self, down_to):
            self.subsampled.add_seqrecord(seq)
        # Clean up #
        self.subsampled.close()
        # Did it work #
        assert len(self.subsampled) == down_to

    def rename_with_num(self, prefix="", new_path=None, remove_desc=True):
        """Rename every sequence based on a prefix and a number"""
        # Temporary path #
        if new_path is None: numbered = self.__class__(new_temp_path())
        else: numbered = self.__class__(new_path)
        # Generator #
        def numbered_iterator():
            for i,read in enumerate(self):
                read.id = prefix + str(i)
                if remove_desc: read.description = ""
                yield read
        # Do it #
        numbered.write(numbered_iterator())
        numbered.close()
        # Replace it #
        if new_path is None:
            os.remove(self.path)
            shutil.move(numbered, self.path)

    def extract_length(self, lower_bound, upper_bound, new_path=None, cls=None):
        """Extract a certain length fraction and place them in a new file"""
        # Temporary path #
        cls = cls or self.__class__
        fraction = cls(new_temp_path()) if new_path is None else cls(new_path)
        # Generator #
        if lower_bound is None: lower_bound = 0
        def fraction_iterator():
            for read in self:
                if lower_bound <= len(read) <= upper_bound:
                    yield read
        # Do it #
        fraction.write(fraction_iterator())
        fraction.close()
        return fraction

    def rename_sequences(self, new_fasta, mapping):
        """Given a new file path, will rename all sequences in the
        current fasta file using the mapping dictionary also provided."""
        assert isinstance(new_fasta, FASTA)
        new_fasta.create()
        for seq in self:
            new_name = mapping[seq.id]
            nucleotides = str(seq.seq)
            new_fasta.add_str(nucleotides, new_name)
        new_fasta.close()

    def extract_sequences(self, new_fasta, ids):
        """Will take all the sequences from the current file who's id appears in
        the ids given and place them in the new file path given."""
        assert isinstance(new_fasta, FASTA)
        new_fasta.create()
        for seq in self:
            if seq.id in ids: new_fasta.add_seq(seq)
        new_fasta.close()

    def align(self, out_path=None):
        """We align the sequences in the fasta file with muscle"""
        if out_path is None: out_path = self.prefix_path + '.aln'
        sh.muscle38("-in", self.path, "-out", out_path)
        return AlignedFASTA(out_path)

    def template_align(self, ref_path):
        # Run it #
        sh.mothur("#align.seqs(candidate=%s, template=%s, search=blast, flip=false, processors=8);" % (self.path, ref_path))
        # Move things #
        shutil.move(self.path[:-6] + '.align', self.aligned_path)
        shutil.move(self.path[:-6] + '.align.report', self.report_path)
        shutil.move(self.path[:-6] + '.flip.accnos', self.accnos_path)
        # Clean up #
        if os.path.exists('formatdb.log'): os.remove('formatdb.log')
        if os.path.exists('error.log') and os.path.getsize('error.log') == 0: os.remove('error.log')
        for p in sh.glob('mothur.*.logfile'): os.remove(p)

    @property_cached
    def graphs(self):
        return dict(((cls,getattr(graphs, cls)(self)) for cls in graphs.__all__))

    @property_cached
    def length_dist(self):
        graph = self.graphs['LengthDist']
        if not graph: graph.plot()
        return graph

    def index(self):
        """Create two indexes. For both bowtie2 and samtools on the contigs fasta file."""
        sh.bowtie2_build(self.contigs_fasta, self.contigs_fasta)
        sh.samtools('faidx', self.contigs_fasta)

    def index_bowtie(self):
        """Create an index on the fasta file compatible with bowtie2"""
        sh.bowtie2_build(self.path, self.path)

    def index_samtools(self):
        """Create an index on the fasta file compatible with samtools"""
        sh.samtools('faidx', self.path)

################################################################################
# Expose objects #
from fasta.fastq import FASTQ
from fasta.aligned import AlignedFASTA
from fasta.paired import PairedFASTQ
from fasta.sizes import SizesFASTA
from fasta.qual import QualFile
from fasta.splitable import SplitableFASTA