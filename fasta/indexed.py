#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Internal modules #
from fasta import FASTA
from plumbing.databases import convert_to_sql
from plumbing.databases.sqlite_database import SQLiteDatabase
from plumbing.common import GenWithLength

# Third party modules #
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tqdm import tqdm

# Constants #
base_keys = ('id', 'description', 'seq')

###############################################################################
class DatabaseFASTA(SQLiteDatabase):

    def __init__(self, path=None):
        self.path = path
        self.factory = lambda cursor, row: SeqRecord(Seq(row[2]), id=row[0], description=row[1])

    def parse(self):
        pass

###############################################################################
def generate_values(path, progress=False):
    seqs = SeqIO.parse(path, 'fasta')
    if not progress:
        for seq in seqs: yield seq.id, seq.description, str(seq.seq)
    if progress:
        for seq in tqdm(GenWithLength(seqs, len(FASTA(path)))):
            yield seq.id, seq.description, str(seq.seq)

###############################################################################
def fasta_to_sql(source, dest):
    values = generate_values(source, progress=True)
    convert_to_sql(dest, base_keys, values)
    return DatabaseFASTA(dest)