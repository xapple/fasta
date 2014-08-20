# Internal modules #
from plumbing.graphs import Graph

# Third party modules #
from matplotlib import pyplot

# Constants #
__all__ = ['LengthDistribution']

################################################################################
class LengthDistribution(Graph):
    """Simple length distribution"""
    short_name = 'length_distribution'
    left = 0.1

    def plot(self, loglog=False):
        # Data #
        counts = self.parent.lengths
        # Plot #
        fig = pyplot.figure()
        pyplot.bar(counts.keys(), counts.values(), 1.0, color='gray', align='center')
        title = 'Distribution of sequence lengths'
        axes = pyplot.gca()
        axes.set_title(title)
        axes.set_xlabel('Length of sequence in nucleotides')
        axes.set_ylabel('Number of sequences with this length')
        axes.yaxis.grid(True)
        # Log log #
        if loglog:
            axes.set_xscale('symlog')
            axes.set_yscale('symlog')
        # Save it #
        self.save_plot(fig, axes, sep=('y'))
        pyplot.close(fig)