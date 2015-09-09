# -*- coding: utf-8 -*-

# Built-in modules #

################################################################################
def add_dummy_scores(iteratable, score=0):
    """Add zero scores to all sequences"""
    for seq in iteratable:
        seq.letter_annotations["phred_quality"] = (score,)*len(seq)
        yield seq