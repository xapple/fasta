#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

################################################################################
def add_dummy_scores(iterable, score=0):
    """Add zero scores to all sequences."""
    for seq in iterable:
        seq.letter_annotations["phred_quality"] = (score,)*len(seq)
        yield seq