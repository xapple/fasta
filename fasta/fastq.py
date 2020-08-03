#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, sys, shutil

# Internal modules #
from fasta import FASTA

# First party modules #
from plumbing.common     import average
from plumbing.cache      import property_cached
from autopaths.tmp_path  import new_temp_path
from autopaths.file_path import FilePath

# Third party modules #
import sh
from Bio import SeqIO

################################################################################
class FASTQ(FASTA):
    """
    A single FASTQ file somewhere in the filesystem.
    It is very similar to a FASTA file object, except we have a few added
    methods and properties for dealing with the quality scores.
    """

    ext    = 'fastq'
    format = 'fastq'

    #----------------------------- Properties --------------------------------#
    @property_cached
    def count(self):
        # Import module #
        from shell_command import shell_output
        # Case when we are not compressed #
        if not self.gzipped:
            return int(int(shell_output("cat %s | wc -l" % self.path)) / 4)
        # If we are gzipped we can just use zcat or gzcat on macOS #
        program = 'gzcat' if sys.platform != 'linux' else 'zcat'
        command = "%s %s | wc -l" % (program, self.path)
        return int(int(shell_output(command)) / 4)

    @property_cached
    def avg_quality(self):
        """
        Returns a single float.
        Quite slow and inefficient computationally.
        Think of a better implementation?
        """
        mean = average(s for r in self for s in r.letter_annotations["phred_quality"])
        self.close()
        return mean

    #----------------------------- Conversion --------------------------------#
    def to_fasta(self, path, verbose=False):
        # Select verbosity #
        import tqdm
        wrapper = tqdm.tqdm if verbose else lambda x: x
        # Do it #
        with open(path, 'w') as handle:
            for r in wrapper(self): SeqIO.write(r, handle, 'fasta')
        # Return #
        return FASTA(path)

    def to_qual(self, path, verbose=False):
        # Select verbosity #
        import tqdm
        wrapper = tqdm.tqdm if verbose else lambda x: x
        # Do it #
        with open(path, 'w') as handle:
            for r in wrapper(self): SeqIO.write(r, handle, 'qual')
        # Return #
        return FilePath(path)

    #-------------------------------- Tools ----------------------------------#
    @property_cached
    def validator(self):
        """Validate the format of this FASTQ."""
        from fasta.validator import Validator
        return Validator(self.path)

    @property_cached
    def fastqc(self):
        """Run the FastQC software on this FASTQ."""
        from fasta.fastqc import FastQC
        return FastQC(self)

    #-------------------------------- PHRED Format ---------------------------#
    def guess_phred_format(self):
        """
        Guess the PHRED score format. The gold standard is the first one, aka
        the one called "Sanger". Sanger encoding is exactly equivalent to
        Illumina-1.8 encoding.
        In other words, they finally gave up with their alternative standards.
        """
        # Possibilities #
        self.intervals = {
            'Sanger':       (33,  74),    # <- This one is the gold standard
            'Solexa':       (59, 105),    # <- This one is rare
            'Illumina-1.3': (64, 105),    # <- These were abandoned after they wised up.
            'Illumina-1.5': (67, 105),    # <- These were abandoned after they wised up.
        }
        # Initialize variables #
        glob_min  = 9999
        glob_max  = 0
        take_next = False
        # Stop after #
        max_sequences_to_consider = 2000000
        # Error message #
        message = "Multiple PHRED encodings possible for '%s'.\n"
        # Loop #
        self.open()
        for i, line in enumerate(self.handle):
            # Are we on the PHRED line? #
            if not take_next:
                take_next = line.startswith('+')
                continue
            # Now we are on the PHRED line! #
            take_next = False
            current_min, current_max = self.get_qual_range(line.rstrip('\n'))
            if current_min < glob_min or current_max > glob_max:
                glob_min = min(current_min, glob_min)
                glob_max = max(current_max, glob_max)
                valid_encodings = [e for e, r in self.intervals.items() if glob_min >= r[0] and glob_max <= r[1]]
                if len(valid_encodings) == 0:
                    message = "Illegal PHRED encoding for '%s'.\n"
                    break
                if len(valid_encodings) == 1:
                    return valid_encodings[0]
            # Stop condition #
            if i >= max_sequences_to_consider: break
        # It didn't work #
        message += "Maximum detected value is: %s.\n"
        message += "Minimum detected value is: %s.\n"
        message += "Possibilities are: %s."
        message  = message % (self, glob_max, glob_min, ' or '.join(valid_encodings))
        raise Exception(message)

    def get_qual_range(self, phred_string):
        """
        >>> self.get_qual_range("DLXYXXRXWYYTPMLUUQWTXTRSXSWMDMTRNDNSMJFJFFRMV")
        (68, 89)
        """
        values = [ord(char) for char in phred_string]
        return min(values), max(values)

    #------------------------- Conversion methods ----------------------------#
    def phred_13_to_18(self, new_path=None, in_place=True):
        """Illumina-1.3 format conversion to Illumina-1.8 format via BioPython."""
        # New file #
        if new_path is None: new_fastq = self.__class__(new_temp_path(suffix=self.extension))
        else:                new_fastq = self.__class__(new_path)
        # Do it #
        self.format = 'fastq-illumina'
        new_fastq.open('w')
        new_fastq.handle.writelines(seq.format('fastq-sanger') for seq in self)
        new_fastq.close()
        self.format = 'fastq-sanger'
        # Return #
        if in_place:
            os.remove(self.path)
            shutil.move(new_fastq, self.path)
            return self
        else: return new_fastq

    def phred_13_to_18_sed(self, new_path=None, in_place=True):
        """Illumina-1.3 format conversion to Illumina-1.8 format via sed (faster)."""
        # String #
        sed_command = r"""4~4y/@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghi/!"#$%&'\''()*+,-.\/0123456789:;<=>?@ABCDEFGHIJ/"""
        # Faster with bash utilities #
        if in_place is True:
            sh.sed('-i', sed_command, self.path)
            return self
        # New file #
        if new_path is None: new_fastq = self.__class__(new_temp_path())
        else:                new_fastq = self.__class__(new_path)
        sh.sed(sed_command + " " + new_fastq, self.path)
        return new_fastq