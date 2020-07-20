#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Internal modules #
from plumbing.graphs import Graph
from autopaths.file_path import FilePath

# Third party modules #
import numpy
from matplotlib import pyplot

# Constants #
__all__ = ['LengthDist', 'LengthHist']

################################################################################
class LengthDist(Graph):
    """The length distribution of the sequences with a bar plot."""

    short_name   = 'length_dist'
    sep          = ('x', 'y')
    y_grid       = True
    width        = 10
    height       = 6
    remove_frame = True

    def __init__(self, parent):
        self.parent = parent
        self.path = FilePath(self.parent.prefix_path + '_len_dist.pdf')

    def plot(self, **kwargs):
        # Data #
        counts = self.parent.lengths_counter
        # Plot #
        fig = pyplot.figure()
        pyplot.bar(counts.keys(), counts.values(), 1.0, color='gray', align='center')
        axes = pyplot.gca()
        # Information #
        title = 'Distribution of sequence lengths'
        axes.set_title(title)
        axes.set_xlabel('Length of sequence in nucleotides')
        axes.set_ylabel('Number of sequences with this length')
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)
        # For convenience #
        return self

################################################################################
class LengthHist(Graph):
    """The length distribution of the sequences with a histogram."""

    short_name   = 'length_hist'
    sep          = ('x', 'y')
    y_grid       = True
    width        = 10
    height       = 6
    remove_frame = True

    def __init__(self, parent):
        self.parent = parent
        self.path = FilePath(self.parent.prefix_path + '_len_hist.pdf')

    def plot(self, bins=80, **kwargs):
        # Data #
        counts = list(self.parent.lengths)
        # Linear bins in logarithmic space #
        if 'log' in kwargs.get('x_scale', ''):
            start, stop = numpy.log10(1), numpy.log10(max(counts))
            bins = list(numpy.logspace(start=start, stop=stop, num=bins))
            bins.insert(0, 0)
        # Plot #
        fig = pyplot.figure()
        pyplot.hist(counts, bins=bins, color='gray')
        axes = pyplot.gca()
        # Information #
        title = 'Histogram of sequence lengths'
        axes.set_title(title)
        axes.set_xlabel('Length of sequence in nucleotides')
        axes.set_ylabel('Number of sequences with this length')
        # X lim #
        axes.set_xlim(min(counts), axes.get_xlim()[1])
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)
        # For convenience #
        return self