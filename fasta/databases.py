# Futures #

# Built-in modules #

# Internal modules #

# Third party modules #
from ftputil import FTPHost

###############################################################################
class RefSeq(object):
    """We need to download the raw sequences of the refseq protein database."""

    ftp_url = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/"

    all_paths = """
    /lorem
    """

    def __init__(self, base_dir=None, seq_type='prot'):
        # Attributes #
        self.seq_type = seq_type
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def download(self):
        """Retrieve all files from the FTP site"""
        self.ftp = FTPHost(self.ftp_url)
        self.ftp_files = self.ftp.listdir()
        files_we_want = files