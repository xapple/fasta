# Built-in modules #

# Internal modules #
from fasta import FASTA

###############################################################################
class QualFile(FASTA):
    """A single QUAL file somewhere in the filesystem"""
    format    = 'qual'
    extension = 'qual'