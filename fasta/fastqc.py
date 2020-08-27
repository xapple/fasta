#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, shutil, multiprocessing

# First party modules #
from fasta import FASTQ
from autopaths.dir_path       import DirectoryPath
from autopaths.tmp_path       import new_temp_dir
from plumbing.cache           import property_cached
from plumbing.check_cmd_found import check_cmd
from plumbing.apt_pkg         import get_apt_packages
from plumbing.scraping        import download_from_url

# Third party modules #
import sh

###############################################################################
class FastQC:
    """
    Takes care of running the FastQC program on a given FASTQ file.
    See http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    Expects version 0.11.9.
    """

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source.path)

    def __init__(self, source, dest=None):
        # Source and destination #
        self.source = FASTQ(source)
        self.dest = DirectoryPath(dest)
        # Default case #
        if dest is None:
            self.dest = DirectoryPath(self.source.prefix_path + '.fastqc')

    #---------------------------- Installing ---------------------------------#
    apt_packages = ['default-jre']
    zip_url = "http://www.bioinformatics.babraham.ac.uk/projects/" \
              "fastqc/fastqc_v0.11.9.zip"

    @classmethod
    def check_installed(cls, exception=True):
        """
        Try to determine if the FastQC software is installed and
        accessible.
        """
        return check_cmd('fastqc', exception, cls.install.__doc__)

    @classmethod
    def install(cls, prefix="~/programs/FastQC/"):
        """
        To automatically download and install the FastQC software on this
        computer and for the current user, type these commands in python:

            >>> from fasta.fastqc import FastQC
            >>> FastQC.install()
        """
        # Start with required apt packages #
        get_apt_packages(cls.apt_packages, verbose=True)
        # Make a temporary directory #
        tmp_dir = new_temp_dir()
        # Download tarball #
        zip_loc = download_from_url(cls.zip_url, tmp_dir, stream=True,
                                    progress=True)
        # Uncompress #
        src_dir = zip_loc.unzip_to(prefix, single=False).sub_directory
        # Set executable permissions #
        bin_loc = src_dir + 'fastqc'
        bin_loc.permissions.make_executable()
        # The directory that contains the executable #
        bin_dir = src_dir.with_tilda[:-1].replace('~', '$HOME')
        # Suggest adding to the $PATH #
        print("\nFastQC was installed successfully. You should now "
              "add this line to your .bash_profile: \n\n    "
              "export PATH=%s:$PATH\n" % bin_dir)

    #---------------------------- Running ------------------------------------#
    def __call__(self, cpus=None):
        # Check it is installed #
        self.check_installed()
        # Check version #
        assert "v0.11.9" in sh.fastqc('-version')
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Make a temporary directory #
        self.tmp_dir = new_temp_dir()
        # Run it #
        sh.fastqc(self.source, '-o', self.tmp_dir, '-t', cpus, '--extract')
        # Get location of results #
        components = self.source.prefix.split('.')
        if components[-1] == 'gz':    components.pop()
        if components[-1] == 'fastq': components.pop()
        # Reassemble the components #
        created_name = '.'.join(components) + '_fastqc/'
        # This will be the name of the directory that fastqc created #
        created_dir  = self.tmp_dir + created_name
        # Move results #
        if self.dest.exists: shutil.rmtree(self.dest)
        shutil.move(created_dir, self.dest)
        self.tmp_dir.remove()
        # Return #
        return self.results

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the FastQC software was run already and the results are
        stored on the filesystem. Return False if it was not yet run.
        """
        return os.path.exists(self.dest + 'Images/per_base_quality.png')

    @property_cached
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from FastQC " \
                  "before running the tool."
            raise Exception(msg)
        # Return the results #
        return FastQCResults(self.dest)

###############################################################################
class FastQCResults(DirectoryPath):
    """A directory with the results from FastQC."""

    all_paths = """
                /Images/per_base_quality.png
                /Images/per_sequence_quality.png
                """

    def __bool__(self): return self.per_base_qual.exists

    @property
    def per_base_qual(self): return self.p.per_base_quality

    @property
    def per_seq_qual(self): return self.p.per_sequence_quality