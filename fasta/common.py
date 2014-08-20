# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from plumbing.autopaths import DirectoryPath, FilePath

###############################################################################
class FastQCResults(DirectoryPath):
    """A directory with the results from FastQC"""

    def __init__(self, directory):
        if not directory.endswith('/'): directory += '/'
        self.path = directory

    @property
    def per_base_qual(self):
        path = FilePath(self.path + 'Images/per_base_quality.png')
        if not path.exists: raise Exception("Attempted to access '%s' which doesn't exist" % path)
        return path

    @property
    def per_seq_qual(self):
        path = FilePath(self.path + 'Images/per_sequence_quality.png')
        if not path.exists: raise Exception("Attempted to access '%s' which doesn't exist" % path)
        return path