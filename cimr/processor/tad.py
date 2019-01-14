
"""checks input file tagged as tad"""

import logging

from cimr.processor.util import findfile

def annotatetad(filename):
    """read the annotation file to integrate into the cimr database
    expects the following information:
    chrom: chromosome id in the format of chr+number
    start: starting position of the tad block
    end: ending position of the tad block
    tadid: unique id of the tad block
    comment: optional column containing additional information
    e.g.
    chr1	1960001	2400000	tad0|hg19|chr1:1960001-2400000	1000"""
    import pandas
    filename = findfile(filename)
    annotdf = pandas.read_csv(filename, sep='\t', header=0)
    return annotdf


