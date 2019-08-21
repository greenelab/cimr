
"""Random walk with restarts for cimr network subprocess"""

import sys
import csv
import numpy
from sklearn.preprocessing import normalize

class Rwalker:
    """Class for multi-graph walk to convergence,
    using matrix computation

    Attributes:
    -----------
    adj_matrix
    restart_prob
    walk_original

    Notes:
    ------
    adopted from https://github.com/greenelab/netwas_rebooted.git

    """

    def __init__(self, delim):
        self.delim = delim

    def _build_matrices(self, original_graph, low_list):
        """
        """
        self.original_graph = original_graph





