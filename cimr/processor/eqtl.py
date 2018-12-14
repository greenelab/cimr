#!/usr/bin/env python3

"""checks input file flagged as eqtl."""

# could delete imports once incorporated into continuous_mr.py
import os
import sys
import csv
import pandas
import pathlib
import logging


# gene : gene id in 
# rsnum : rs id of the variant
# constant_id : chromosome\_position\_referenceallele\_alternateallele\_genomebuild
# e.g. chr2_128747549_G_T_b38
# inc_allele : allele with respect to which variant's effect sizes are estimated
# inc_afrq : allele frequency of inc_allele
# beta : beta coefficient estimate for the association effect of the variant on the gene 
# se : standard error of the beta
# pval : p-value of the beta estimate


def findfile(filename):
    """find the file indicated in the prompt"""
    if(pathlib.Path(filename).resolve()):
        filename = str(filename)
        return sumdata
    else:
        logging.error(f'no file {filename} found for processing.')
        exit()


def getpos(sumdata):
    """check variant_id column and make chrom, pos, ref, alt, build columns"""
    temp = sumdata.variant_id.str.split('_', expand=True)
    if not temp.empty:
        sumdata['chrom'] = temp[0]
        sumdata['pos'] = temp[1]
        sumdata['ref'] = temp[2]
        sumdata['alt'] = temp[3]
        sumdata['build'] = temp[4]
    return sumdata


def check_chrom(sumdata):
    """assumes numeric coding but will allow chr+number
    - check for autosomal chromosomes
    - change if different from the specified format
    - discard non-autosomal chromosomes from main input
    """
    chromdict = {str(i):'chr'+str(i) for i in range(1,23)}
    chromstr = ['chr'+str(i) for i in range(1,23)]
    chromint = [i for i in range(1,23)]
    chroms = sumdata.chrom.drop_duplicates().values
    if len(chroms) > 21:
        pass
    else:
        logging.warning(f'input file does not include all autosomal chromosomes.')
        logging.warning(f'chromosome(s) included in the input file: %s'%(chroms,))
    if len(set(chroms) & set(chromstr)) > 1:
        sumdata = sumdata[sumdata['chrom'].isin(chromstr)]
    elif len(set(chroms) & set(chromint)) > 1:
        sumdata['chrom'] = sumdata['chrom'].map(chromdict)
        sumdata = sumdata[sumdata['chrom'].isin(chromstr)]
    else:
        logging.error(f'chromosome id needs to be checked.')
    remainder = list(set(chroms) - set(chromstr) - set(chromint))
    logging.warning(f'chromosome(s) not used for analysis using the current version of cimr: %s'%(remainder,))
    return sumdata


def check_numeric(sumdata):
    """beta, se, pval columns are expected to be numeric"""
    from pandas.api.types import is_numeric_dtype
    for col in ['beta', 'se', 'pval']:
        try:
            if is_numeric_dtype(sumdata[col]):
                print(col+' is numeric')
                return sumdata
            else:
                numdata = (sumdata
                           .drop([col], axis=1)
                           .join(sumdata[col].apply(pandas.to_numeric, errors='coerce')))
                numcol = numdata[col].isnull().values().sum()
                logging.error(f'%s rows in %s are non-numeric'%(numcol,col,))
                return numdata
        except:
            logging.error(f'the format of %s is not testable.'%(col,))
            exit()

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


def readfile(filename):
    """read the input file as a pandas dataframe. check if empty"""
    filename = findfile(filename)
    sumdata = pandas.read_csv(filename, sep='\t', header=0)
    if not sumdata.empty:
        sumdata = check_chrom(getpos(sumdata))
    else:
        logging.error(f'no content in uploaded file {filename}.')
    return sumdata


if __name__ == '__main__':
    filename = sys.argv[1]
    sumdata = readfile(filename)


    
