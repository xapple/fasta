# Futures #

# Built-in modules #
import fnmatch

# Internal modules #
from plumbing.autopaths import AutoPaths

# Third party modules #
from ftputil import FTPHost
from tqdm import tqdm

# Constants #
directory = "/glob/lucass/"

###############################################################################
class RefSeqBacteriaProtNR(object):
    """the RefSeq sequences for only bacteria, only protein, and only the
    non-redundant version. We will download the raw sequences by FTP
    and format them as a blast data base."""

    short_name = "refseq_bact_prot_nr"
    ftp_url    = "ftp.ncbi.nlm.nih.gov"
    ftp_dir    = "/refseq/release/bacteria/"

    all_paths = """
    /raw/
    /blast_db/
    """

    def __init__(self):
        # Attributes #
        self.seq_type = 'prot'
        self.base_dir = directory + self.short_name
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def download(self):
        """Retrieve all files from the FTP site"""
        self.ftp = FTPHost(self.ftp_url, "anonymous")
        self.ftp.chdir(self.ftp_dir)
        self.ftp_files = self.ftp.listdir(self.ftp.curdir)
        self.pattern = 'bacteria.nonredundant_protein.*.protein.faa.gz'
        self.files_we_want = [f for f in self.ftp_files if fnmatch.fnmatch(f, self.pattern)]
        for f in tqdm(self.files_we_want): self.ftp.download(f, self.p.raw_dir + f)

###############################################################################
refseq_bact_prot_nr = RefSeqBacteriaProtNR()