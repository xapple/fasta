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
    sep        = ('x')
    x_grid     = True

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