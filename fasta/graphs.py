# Internal modules #
from plumbing.graphs import Graph
from plumbing.autopaths import FilePath

# Third party modules #
from matplotlib import pyplot

# Constants #
__all__ = ['LengthDist']

################################################################################
class LengthDist(Graph):
    """The length distribution of the sequences"""
    short_name = 'length_dist'

    def __init__(self, parent):
        self.parent = parent
        self.path = FilePath(self.parent.prefix_path + '_len_dist.pdf')

    def plot(self, x_log=False, y_log=False):
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
        axes.xaxis.grid(False)
        # Add logarithm to axes #
        if x_log: axes.set_xscale('symlog')
        if y_log: axes.set_yscale('symlog')
        # Save it #
        self.save_plot(fig, axes, sep=('x'))
        # For convenience #
        return self