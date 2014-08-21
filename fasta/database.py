# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from fasta import FASTA
from plumbing.database import Database, convert_to_sql

# Third party modules #
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Constants #
base_keys = ('id', 'description', 'seq')

###############################################################################
class DatabaseFASTA(Database, FASTA):

    def __init__(self, path=None):
        self.path = path
        self._factory = None
        #self.factory = lambda cursor, row: SeqRecord(Seq(row[2]), id=row[0], description=row[1])

    def parse(self):
        self.open()
        return SeqIO.parse(self.handle, self.format)

###############################################################################
def fasta_to_sql(source, dest):
    def generate_values():
        seqs = SeqIO.parse(source, 'fasta')
        for seq in seqs: yield (seq.id, seq.description, str(seq.seq))
    convert_to_sql(source, dest, base_keys, generate_values())
    return DatabaseFASTA(dest)