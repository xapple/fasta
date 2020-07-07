#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import sys

# First party modules #
from autopaths.tmp_path import new_temp_dir
from plumbing.apt_pkg   import get_apt_packages, check_apt_exists
from plumbing.scraping  import download_from_url
from plumbing.check_cmd_found import check_cmd

# Third party modules #
import sh

###############################################################################
class Validator:
    """
    Determine if a FASTQ file is valid by calling
    https://github.com/statgen/fastQValidator
    on the file.
    Documentation is at:
    https://genome.sph.umich.edu/wiki/FastQValidator
    """

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.path)

    def __init__(self, path):
        self.path = path

    #---------------------------- Installing ---------------------------------#
    apt_packages = ['g++', 'libssl-dev', 'zlib1g-dev']

    @classmethod
    def check_installed(cls, exception=True):
        """
        Try to determine if the fastQValidator software is installed and
        accessible on this machine.
        """
        return check_cmd('fastQValidator', exception, cls.install.__doc__)

    @classmethod
    def install(cls, prefix="~/programs/fastQValidator/"):
        """
        To automatically download and install the fastQValidator software
        on this computer and for the current user, type these commands in
        python:

            >>> from fasta.validator import Validator
            >>> Validator.install()
        """
        # Check we are on an OS with aptitude #
        check_apt_exists()
        # Start with the required apt packages #
        get_apt_packages(cls.apt_packages, verbose=True)
        # Download tarball 1 #
        tmp_dir_1 = new_temp_dir()
        tgz_url_1 = 'https://github.com/statgen/libStatGen/archive/master.tar.gz'
        tgz_loc_1 = download_from_url(tgz_url_1, tmp_dir_1, stream=True, progress=True)
        src_dir_1 = tgz_loc_1.untargz_to()
        # Download tarball 2 #
        tmp_dir_2 = new_temp_dir()
        tgz_url_2 = 'https://github.com/statgen/fastQValidator/archive/master.tar.gz'
        tgz_loc_2 = download_from_url(tgz_url_2, tmp_dir_2, stream=True, progress=True)
        src_dir_2 = tgz_loc_2.untargz_to()
        # Uncompressed 1 #
        src_dir_1 = src_dir_1.sub_directory
        # Uncompressed 2 #
        src_dir_2 = src_dir_2.sub_directory
        # Make 1 #
        sh.make('-C', src_dir_1, _out=sys.stdout, _err=sys.stderr)
        # Make 2 #
        sh.make('-C', src_dir_2, 'LIB_PATH_FASTQ_VALIDATOR=%s' % src_dir_1,
                  _out=sys.stdout, _err=sys.stderr,)
        # Move the executable #
        binary = src_dir_2 + 'bin/fastQValidator'
        path = binary.move_to(prefix, overwrite=True)
        # The directory that contains the executable #
        bin_dir = path.directory.with_tilda[:-1].replace('~', '$HOME')
        # Suggest adding to the $PATH #
        print("\nfastQValidator was installed successfully. You should now "
              "add this line to your .bash_profile: \n\n    "
              "export PATH=%s:$PATH\n" % bin_dir)

    #---------------------------- Running ------------------------------------#
    def __call__(self, exception=True):
        # Check it is installed #
        self.check_installed()
        # Run software #
        result = sh.fastQValidator('--file', self.path)
        # Check result #
        if "FASTQ_SUCCESS" not in result:
            msg = "The fastq file '%s' failed to validate."
            msg = msg % self.path
            if exception: raise Exception(msg)
            return False
        # Return #
        return True