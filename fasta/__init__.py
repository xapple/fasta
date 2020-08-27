#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, sys

# Constants #
url         = 'http://xapple.github.io/fasta/'
repo_url    = 'http://github.com/xapple/fasta/'
__version__ = '2.1.8'

# Get paths to module #
self       = sys.modules[__name__]
module_dir = os.path.dirname(self.__file__)

# Expose objects for convenience #
from fasta.core      import FASTA
from fasta.fastq     import FASTQ
from fasta.aligned   import AlignedFASTA
from fasta.paired    import PairedFASTQ, PairedFASTA
from fasta.sizes     import SizesFASTA
from fasta.qual      import QualFile
from fasta.splitable import SplitableFASTA