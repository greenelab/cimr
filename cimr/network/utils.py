
"""Utility functions for cimr network subprocess

(c) YoSon Park
"""

import sys
import numpy
import logging


class EdgeNetworker:
    """Utilities for edge list input files

    Parameters:
    -----------

    Notes:
    ------

    """
    def __init__(self, filename, edgetype='binary', diags='1.0'):
        self.filename = filename
        self.edgetype = edgetype
        self.diags = diags


    def transform_edgelist(self):
        """Reads a list of edges and returns a numpy binary
        Parameters:
        -----------
        node1\tnode2\t(weight)

        Notes:

        """
        self.data = []
        self.nodemap = {}
        current_node = 0

        try:
            with open(self.filename, 'r') as f:
                for line in f:
                    if len(line.split('\t')) == 3:
                        g1, g2, weight = line.split('\t')
                        weight = float(weight)
                        self.edgetype = 'weighted'
                    elif len(line.split('\t')) == 2:
                        g1, g2 = line.split('\t')
                        weight = 1.0
                        self.edgetype = 'binary'
                    else:
                        logging.error(f' not correct number of columns in the file.')
                        sys.exit()

                    if g1 not in self.nodemap:
                        self.nodemap[g1] = current_node
                        # to accomodate liblinear, index starts with 1
                        self.data.append({self.nodemap[g1]: 1})
                        current_node += 1
                    if g2 not in self.nodemap:
                        self.nodemap[g2] = current_node
                        self.data.append({self.nodemap[g2]: 1})
                        current_node += 1

                    self.data[self.nodemap[g1]][self.nodemap[g2]] = weight
                    self.data[self.nodemap[g2]][self.nodemap[g1]] = weight

            self.nodes = sorted(self.nodemap, key=self.nodemap.get)

            return self

        except:
            logging.error(f' the file {self.filename} is not readable.')
            sys.exit()


    def make_adj_matrix(self, outfile):
        """Make numpy array from the edgelist"""
        self.outfile = outfile
        matrixscale = len(self.nodes)
        self.adjmatrix = numpy.zeros(matrixscale, matrixscale)

        for idx, _ in enumerate(self.data):
            col_inds = list(self.data[idx].keys())
            col_weights = list(self.data[idx].values())
            self.adjmatrix[idx,col_inds] = col_weights
        numpy.fill_diagonal(self.adjmatrix, self.diags)
        numpy.save(self.outfile+'.npy', self.adjmatrix)
        numpy.savetxt(self.outfile+'_nodelist.txt', self.nodes, fmt='%s')

        return self


class BaseNetworker:
    """Utilities for matrix input files used in network analysis

    Parameters:
    -----------

    Notes:
    ------
    """
    def __init__(self, filename):
        self.filename = filename+'.npy'
        self.nodefile = filename+'_nodelist.txt'


    def load_adjmatrix(self, diags=False):
        self.adjmatrix = numpy.load(self.filename, delimiter='\t')

        if diags:
            numpy.fill_diagonal(self.adjmatrix, 1)

        self.nodes = numpy.loadtxt(self.nodefile)
