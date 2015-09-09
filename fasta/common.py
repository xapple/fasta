# -*- coding: utf-8 -*-

# Built-in modules #

################################################################################
def add_dummy_scores(iteratable):
    """Add zero scores to all sequences"""
    for seq in iteratable:
        seq.letter_annotations["phred_quality"] = (0,)*len(seq)
        yield seq