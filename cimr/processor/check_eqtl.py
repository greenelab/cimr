#!/usr/bin/env python3

"""checks input file from contributors flagged as eqtl."""

# could delete imports once incorporated into continuous_mr.py
import os
import sys
import pandas

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


def readfile(infile):
    pass

def check_chrom(infile):
    pass

