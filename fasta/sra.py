#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Built-in modules #
import sys, os, functools, uuid, inspect, multiprocessing

# First party modules #
from autopaths import Path
from autopaths.file_path import FilePath
from autopaths.tmp_path import new_temp_dir
from plumbing.check_cmd_found import check_cmd
from plumbing.scraping import download_from_url
from fasta import FASTQ, PairedFASTQ

# Third party modules #
import sh

###############################################################################
class DumpSRA:
    """
    Takes care of running the `sra-toolkit` program to extract a FASTQ (or a
    pair of FASTQs) from an SRA file. See:

    https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

    Until Ubuntu 19, `sra-toolkit` was a package in the apt universe. This
    is not the case for Ubuntu 20. See:

    https://askubuntu.com/questions/1232028/sra-toolkit-for-ubuntu-20-04-tls

    If you are on macOS you can just type "brew install sratoolkit".
    """

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, dest=None):
        # Source is the SRA on which the sra-toolkit will be run #
        self.source = FilePath(source)
        # Destination is a FASTQ file that will contain the results #
        self.dest = dest
        # Or alternatively it is a FASTQ pair of files (forward, reverse) #
        if self.dest is None:
            self.dest = self.source.prefix_path + '.fastq'
        # Check it is a FASTQ type class #
        if not isinstance(self.dest, (FASTQ, PairedFASTQ)):
            self.dest = FASTQ(dest)

    #---------------------------- Installing ---------------------------------#
    url = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/" \
          "sratoolkit.3.0.0-ubuntu64.tar.gz"

    @classmethod
    def check_installed(cls, exception=True):
        """
        Try to determine if the `sra-toolkit` software is installed and
        accessible.
        """
        msg = inspect.getdoc(cls.install)
        return check_cmd('fasterq-dump', exception, msg)

    @classmethod
    def install(cls, install_dir="~/programs/sra-toolkit/"):
        """
        To automatically download and install the `sra-toolkit` software on
        this computer and for the current user, type these commands in python:

            >>> from fasta.sra import DumpSRA
            >>> DumpSRA.install()
        """
        # Message #
        print("Installing sra-toolkit into '%s'." % install_dir)
        # Make a temporary directory #
        tmp_dir = new_temp_dir(prefix='sra_install-')
        # Download tarball #
        zip_loc = download_from_url(cls.url, tmp_dir, stream=True,
                                    progress=True)
        # Uncompress #
        zip_loc.untargz_to(tmp_dir)
        src_dir = tmp_dir.sub_directory
        src_dir.move_to(install_dir)
        # The directory that contains the executable #
        bin_dir = src_dir + 'bin/'
        # Mandatory configuration #
        cls.vdb_config_workaround()
        # Suggest adding to the $PATH #
        print("\nThe sra-toolkit was installed successfully. You should now "
              "add this line to your `.bash_profile`: \n\n    "
              "export PATH=%s:$PATH\n" % bin_dir)
        # Add it to the path #
        os.environ["PATH"] += os.pathsep + bin_dir

    @classmethod
    def vdb_config_workaround(cls):
        """
        This method needs to be run to avoid an error when running any of
        the `sra-toolkit` tools. This is due to some very strange and
        nonsensical change introduced in version 2.10.3
        See: https://github.com/ncbi/sra-tools/issues/291
        """
        # The directory where the configuration file must be stored #
        ncbi_settings = Path("~/.ncbi/user-settings.mkfg")
        ncbi_settings.directory.create_if_not_exists()
        # Generate a new random ID #
        new_guid = uuid.uuid4()
        # Add this to the file #
        ncbi_settings.write('/LIBS/GUID = "%s"\n' % new_guid)

    #------------------------------ Running ----------------------------------#
    def __call__(self, cpus=None, verbose=True, progress=True):
        # Message #
        if verbose: print("Running `fasterq-dump` on '%s'" % self.source)
        # Check it is installed #
        self.check_installed()
        # Some questionable change they did in recent versions #
        self.vdb_config_workaround()
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Get the command #
        fasterq = sh.Command("fasterq-dump")
        # Make a temporary directory #
        tmp_dir = new_temp_dir(prefix='sra_dump-')
        # Command line options #
        options = {'O': tmp_dir,  # --outdir
                   't': tmp_dir,  # --temp
                   'e': cpus}     # --threads
        # Extras #
        if progress:
            print("There are two steps. First 'join', then 'concat'.")
            options['progress'] = True
            options['_out'] = sys.stdout
            options['_err'] = sys.stderr
        # Run it #
        fasterq(self.source, **options)
        # Is it a single FASTQ we are expecting or a pair? #
        if isinstance(self.dest, PairedFASTQ):
            fwd = tmp_dir + self.source.prefix + '_1.fastq'
            rev = tmp_dir + self.source.prefix + '_2.fastq'
            fwd.move_to(self.dest.rev, overwrite=True)
            rev.move_to(self.dest.fwd, overwrite=True)
        else:
            result = tmp_dir + self.source.prefix + '.fastq'
            result.move_to(self.dest, overwrite=True)
        # Remove the temporary directory #
        tmp_dir.remove()
        # Do we want to compress the result #
        if self.dest.gzipped:
            self.dest.compress(new_path=False, verbose=verbose)
        # Return #
        return self.dest

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the `sra-toolkit` software was run already and the
        results are stored on the filesystem. Return False if it was not yet
        run.
        """
        return self.dest.exists

    @functools.cached_property
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from `sra-toolkit` " \
                  "before running the tool."
            raise Exception(msg)
        # Return the results #
        return self.dest