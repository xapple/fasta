# Built-in modules #
import fnmatch
from collections import OrderedDict

# Internal modules #
from plumbing.autopaths import AutoPaths, FilePath
from plumbing.cache import property_cached

# Third party modules #
from ftputil import FTPHost
from tqdm import tqdm

# Constants #
directory = "/glob/lucass/"

###############################################################################
class Database(object):
    """General database object to inherit from."""
    all_paths = """
    /raw/
    /blast_db/
    """

    def __init__(self, seq_type):
        # Attributes #
        self.seq_type = seq_type
        self.base_dir = directory + self.short_name
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def ftp(self):
        """If the data is to be obtained by FTP here is the ftputil object."""
        ftp = FTPHost(self.ftp_url, "anonymous")
        ftp.chdir(self.ftp_dir)
        return ftp

    @property_cached
    def files_to_retrive(self):
        """The files we want to download with their destinations."""
        files = self.ftp.listdir(self.ftp.curdir)
        return OrderedDict((f, FilePath(self.p.raw_dir+f)) for f in files
                            if fnmatch.fnmatch(f, self.pattern))

    @property_cached
    def files_remaining(self):
        """The files we haven't downloaded yet based on size checks."""
        return OrderedDict((source,dest) for source,dest in self.files_to_retrive.items()
                           if dest.count_bytes != self.ftp.path.getsize(source))

    @property
    def raw_files(self):
        """The files we have downloaded."""
        return self.p.raw.contents

    def download(self):
        """Retrieve all files from the FTP site"""
        for source,dest in tqdm(self.files_remaining.items()):
            dest.remove()
            self.ftp.download(source, dest)

###############################################################################
class RefSeqBacteriaProtNR(Database):
    """the RefSeq sequences for only bacteria, only protein, and only the
    non-redundant version. We will download the raw sequences by FTP
    and format them as a blast data base."""

    short_name = "refseq_bact_prot_nr"
    ftp_url    = "ftp.ncbi.nlm.nih.gov"
    ftp_dir    = "/refseq/release/bacteria/"
    pattern    = 'bacteria.nonredundant_protein.*.protein.faa.gz'

###############################################################################
class RefSeqArchaeaProtNR(Database):
    """the RefSeq sequences for only bacteria, only protein, and only the
    non-redundant version. We will download the raw sequences by FTP
    and format them as a blast data base."""

    short_name = "refseq_arch_prot_nr"
    ftp_url    = "ftp.ncbi.nlm.nih.gov"
    ftp_dir    = "/refseq/release/archaea/"
    pattern    = 'archaea.nonredundant_protein.*.protein.faa.gz'

###############################################################################
refseq_bact_prot_nr = RefSeqBacteriaProtNR('prot')
refseq_arch_prot_nr = RefSeqArchaeaProtNR('prot')