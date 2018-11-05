#!/usr/bin/env python3

"""checks input file flagged as eqtl."""

# could delete imports once incorporated into continuous_mr.py
import os
import sys
import pandas
import pathlib
import logging

# note to self

# gene : gene id in 
# rsnum : rs id of the variant
# constant_id : chromosome\_position\_referenceallele\_alternateallele\_genomebuild
# e.g. chr2_128747549_G_T_hg19
# inc_allele : allele with respect to which variant's effect sizes are estimated
# inc_afrq : allele frequency of inc_allele
# beta : beta coefficient estimate for the association effect of the variant on the gene 
# se : standard error of the beta
# pval : p-value of the beta estimate

def findfile(filename):
    if(pathlib.Path(filename).resolve()):
        filename = str(filename)
        return filename
    else:
        logging.error(f'no file {filename} found for processing.\n')

def readfile(filename):
    filename = findfile(filename)
    sumdata = pandas.read_csv(filename, sep='\t', header=0)
    if not sumdata.empty:
        return sumdata
    else:
        logging.error(f'no content in uploaded file {filename}.\n')
    return sumdata

def check_chrom(infile):
    pass


if __name__ == '__main__':
    filename = sys.argv[1]
    filename = findfile(filename)
    sumdata = readfile(filename)
    print(filename)
    print(sumdata)

    
