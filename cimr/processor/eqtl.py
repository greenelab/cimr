#!/usr/bin/env python3

"""checks input file flagged as eqtl."""

# could delete imports once incorporated into continuous_mr.py
import os
import sys
import csv
import pandas
import pathlib
import logging


def find_file(filename):
    """check if file exists"""
    if(pathlib.Path(filename).exists()):
        filename = str(pathlib.Path(filename).resolve(strict=True))
        return filename
    else:
        logging.error(f'no file {filename} found for processing.\n')

def split_constantid(sumdata):
    """parse the variant_id column for variant information
    variant_id column assumes the following format:
    chromosome_position_referenceallele_alternateallele_genomebuild
    e.g. chr2_128747549_G_T_hg19
    """
    splitdata = sumdata.variant_id.str.split('_').tolist()
    parsedvar = pandas.DataFrame(splitdata, columns=['chrom', 'pos_hg38', 'ref', 'alt', 'build'])
    sumdata = sumdata.join(parsedvar)
    return sumdata

def check_header(sumdata):
    """check header of the input file then return matches
    following are required:
    - gene_id : gene id in 
    - rsnum : rs id of the variant
    - variant_id : variant id in the following format
        chromosome_position_referenceallele_alternateallele_genomebuild
    - pval : p-value of the beta estimate
    following are recommended:
    - beta : beta coefficient estimate for the association effect of the variant on the gene 
    - se : standard error of the beta
    - zscore: z-score if no beta/se information is present (e.g. imputed summary statistic)
    following are optional:
    - tss_distance : distance to the transcription start site of the gene_id
    - ma_samples : samples with minor alleles
    - maf : minor allele frequency
    - inc_allele : allele with respect to which variant's effect sizes are estimated
    - inc_afrq : allele frequency of inc_allele
    following are given by parsing variant_id column:
    - chrom : chromosome id
    - pos_hg38 : genomic position in GRCh38
    - ref : reference allele
    - alt : alternate allele
    - build : genomic build version number for the variant_id
    """
    maxhead = ['gene_id', 'rsnum', 'variant_id', 'pval', 'beta', 'se', 'zscore', 'tss_distance', 'ma_samples', 'maf', 'chrom', 'pos_hg38', 'ref', 'alt', 'build', 'inc_allele', 'inc_afrq']
    sumhead = sumdata.columns
    inchead = list(set(maxhead) & set(sumhead))
    
    if 'variant_id' in inchead:
        pass
    else:
        logging.error(f'no variant information provided.')
        sys.exit()
    return inchead

def check_chrom(chrom):
    """check column containing chromosome IDs based on chrom column in the parsed variant_id column.
    - currently, cimr does not accept sex (23, X, Y) or mitochondrial chromosomes (24, MT, M)
    """
    chromlist = chrom.unique()
    chromcount = len(chromlist)
    logging.info(f'there are {chromcount} unique chromosomes represented in the input file.')
    return chromlist

def read_file(filename):
    """reading as a pandas dataframe
    assumes tab-delimited file with a header"""
    filename = find_file(filename)
    sumdata = pandas.read_csv(filename, sep='\t', header=0)
    if not sumdata.empty:
        sumdata = split_constantid(sumdata)
        chromlist = check_chrom(sumdata.chrom)
        return sumdata
    else:
        logging.error(f'no content in uploaded file {filename}.\n')
    # return sumdata

if __name__ == '__main__':
    filename = sys.argv[1]
    sumdata = readfile(filename)
    sumdata = split_constantid(sumdata)
    check_chrom(sumdata.chrom)

    
