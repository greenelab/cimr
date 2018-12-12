#!/usr/bin/env python

# YoSon Park

""" visualize similarities between networks
1. find random subset of network edges
2. subset all cell type specific networks for 1
3. cosine similarities between cell types
"""


import os
import sys
import csv
import argparse


def select_random_edges(celltype, filesize, offsets):
    """select random edges from the celltype dat file
    assumes a tab-delimited file
    """
    import random
    net = open(celltype+'.dat')
    offset = random.randrange(filesize)
    net.seek(offset)
    # net.readline()
    random_edge = net.readline()
    if len(random_edge) == 0:
        net.seek(0)
        random_edge = net.readline()
    random_edge = random_edge.rstrip().split('\t')
    return offset, random_edge

def write_subset(celltype, filesize, maxcount):
    """write subset of edges chosen by select_random_edges"""
    offsets = []
    edges = []
    while len(edges) < maxcount:
        offset, edge = select_random_edges(celltype, filesize, offsets)
        if not edge in edges:
            offsets.append(offset)
            edges.append(edge)
    with open(celltype+'_'+str(maxcount)+'edges.tsv', 'w') as out:
        csvwriter = csv.writer(out,delimiter='\t')
        csvwriter.writerows(edges)
    return 0

def select_known_edges(celltype, subsetfile):
    """select edges based on a list of edges in the subset file"""
    import pandas
    netfile = celltype+'.dat.gz'
    suffix = subsetfile.replace('global','')
    edgesfile = celltype+suffix
    if os.path.exists(netfile):
        netdf = pandas.read_csv(netfile, sep='\t', compression='gzip', header=None, names=['node0', 'node1', 'edge'], dtype=object)
    if os.path.exists(subsetfile):
        nodedf = pandas.read_csv(subsetfile, sep='\t', header=None, names=['node0', 'node1', 'edge_old'], dtype=object)
    if not netdf.empty and not nodedf.empty:
        edges = netdf.merge(nodedf, on=['node0', 'node1'], left_index=False, right_index=False, how='inner')
        edges = edges[['node0', 'node1', 'edge']]
        edges.to_csv(edgesfile, sep='\t', header=False, index=False)
    return 0

def cosinesim(cellfile, suffix):
    """given a file containing a list of celltypes, 
    estimate cosine similarities between celltypes and save file
    assumes that the file names are celltype + suffix
    resulting file can be used for a heatmap"""
    from scipy import spatial
    import pandas
    # import seaborn
    import itertools
    celltypelist = pandas.read_csv(cellfile, sep='\t', header=0, names=['celltype', 'network'])
    celltypelist = celltypelist.celltype.tolist()
    cellmat = pandas.DataFrame(index=celltypelist, columns=celltypelist)
    for cell0, cell1 in list(itertools.product(celltypelist, celltypelist)):
        cell0edge = pandas.read_csv(cell0+suffix, sep='\t', header=None, names=['node0', 'node1', 'edge'])
        cell1edge = pandas.read_csv(cell1+suffix, sep='\t', header=None, names=['node0', 'node1', 'edge'])
        cellsim = 1 - spatial.distance.cosine(cell0edge.edge, cell1edge.edge)
        cellmat.loc[cell0,cell1] = cellsim
        cellmat.loc[cell1,cell0] = cellsim
        print(cell0, cell1)
    outfile = cellfile.split('.')[0] + '_cossim.tsv'
    # pngfile = outfile.replace('.tsv','.png')
    cellmat.to_csv(outfile, sep='\t', index=True, header=True)
    # seaborn.clustermap(cellmat, cmap='Blues').savefig(pngfile)
    return 0

def parse_arguments():
    """basic arguments
    includes options to change:
    cell type of network (assumed to be the file name)
    filesize to set the randrange stopping point
    maxcount of edges to subset
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--celltype',
        dest='celltype',
        default='global',
        nargs=1,
        type=str,
        help='cell or tissue type that both represents the specificity '
             'of the given network and the file name'
             'e.g. for giant2, file name is assumed to be celltype.dat',
    )
    parser.add_argument(
        '--filesize',
        dest='filesize',
        default=10000000,
        nargs=1,
        type=int,
        help='size of the file containing the network in the format of '
             'edge0 edge1 weight',
    )
    parser.add_argument(
        '--maxcount',
        dest='maxcount',
        default=100000,
        nargs=1,
        type=int,
        help='count of edges to subset from the total network',
    )
    parser.add_argument(
        '--select',
        default=False,
        action='store_true',
        help='select edges based on a list of nodes'
    )
    parser.add_argument(
        '--subsetfile',
        dest='subsetfile',
        default=None,
        nargs=1,
        type=str,
        help='name of the file containing a list of nodes to select '
             'required if --select argument is used.'
    )
    parser.add_argument(
        '--cossim',
        dest='cossim',
        default=False,
        action='store_true',
        help='estimate cosine similarities for a heatmap'
    )
    parser.add_argument(
        '--cellfile',
        dest='cellfile',
        default=None,
        nargs=1,
        type=str,
        help='name of the file contaning a list of celltypes to '
             'estimate cosine similarities for a heatmap'
    )
    args = parser.parse_args()
    return args



if __name__ == '__main__':
    args = parse_arguments()
    celltype = args.celltype[0]
    filesize = args.filesize
    maxcount = args.maxcount
    if args.select:
        if args.subsetfile is not None:
            subsetfile = args.subsetfile[0]
            select_known_edges(celltype, subsetfile)
    if args.cossim:
        if args.cellfile is not None:
            cosinesim(args.cellfile[0], '_100000edges.tsv')
    else:
        write_subset(celltype, filesize, maxcount)




