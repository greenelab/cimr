#!/usr/bin/env python3

"""utilities and common file checks used across different functions in
processor subunits"""

import sys
import pandas
import pathlib
import logging

def findfile(filename):
    """find the file indicated in the prompt"""
    if(pathlib.Path(filename).resolve()):
        filename = str(filename)
        return filename
    else:
        logging.error(f'no file {filename} found for processing.')
        exit()


def getpos(summary_data):
    """check variant_id column and make chrom, pos, ref, alt, build columns"""
    temp = summary_data.variant_id.str.split('_', expand=True)
    if not temp.empty:
        summary_data['chrom'] = temp[0]
        summary_data['pos'] = temp[1]
        summary_data['ref'] = temp[2]
        summary_data['alt'] = temp[3]
        summary_data['build'] = temp[4]
    return summary_data


def checkchrom(summary_data):
    """assumes numeric coding but will allow chr+number
    - check for autosomal chromosomes
    - change if different from the specified format
    - discard non-autosomal chromosomes from main input
    """
    chromdict = {str(i):'chr'+str(i) for i in range(1,23)}
    chromstr = ['chr'+str(i) for i in range(1,23)]
    chromint = [i for i in range(1,23)]
    chroms = summary_data.chrom.drop_duplicates().values
    if len(chroms) > 21 and len(chroms) < 27:
        pass
    elif len(chroms) <= 21:
        logging.warning(f'input file does not include all autosomal chromosomes.')
        logging.warning(f'chromosome(s) included in the input file: %s'%(chroms,))
    else:
        logging.warning(f'input file more than 22 chromosomes.')
        logging.warning(f'chromosome(s) included in the input file: %s'%(chroms,))    
    if len(set(chroms) & set(chromstr)) > 1:
        summary_data = summary_data[summary_data['chrom'].isin(chromstr)]
    elif len(set(chroms) & set(chromint)) > 1:
        summary_data['chrom'] = summary_data['chrom'].map(chromdict)
        summary_data = summary_data[summary_data['chrom'].isin(chromstr)]
    else:
        logging.error(f'chromosome id needs to be checked.')
    remainder = list(set(chroms) - set(chromstr) - set(chromint))
    if len(remainder) > 0:
        logging.warning(f'chromosome(s) not used for analysis using the current version of cimr: %s'%(remainder,))
    return summary_data


def checkhead(summary_header):
    """check header of the input file then return matches
    following are required:
    - gene_id : gene id. expected for eqtls
    - rsnum : rs id of the variant
    - variant_id : variant id in the following format
        chromosome_position_referenceallele_alternateallele_genomebuild
    - pval / pvalue: p-value of the beta estimate
    following are recommended:
    - beta / effect_size: beta coefficient estimate for the association effect of the variant on the gene 
    - se / standard_error: standard error of the beta
    - zscore: z-score if no beta/se information is present (e.g. imputed summary statistic)
    following are optional:
    - tss_distance : distance to the transcription start site of the gene_id
    - ma_samples : samples with minor alleles
    - maf : minor allele frequency
    - inc_allele : allele with respect to which variant's effect sizes are estimated
    - inc_afrq : allele frequency of inc_allele
    following are given by parsing variant_id column:
    - chrom : chromosome id
    - pos: genomic position
    - ref : reference allele
    - alt : alternate allele
    - build : genomic build version number for the variant_id
    """
    maxhead = ['gene_id', 'rsnum', 'variant_id', 'pvalue', 'effect_size', 'standard_error', 'zscore', 'tss_distance', 'ma_samples', 'maf', 'chrom', 'pos', 'ref', 'alt', 'build', 'inc_allele', 'inc_afrq']
    included_header = list(set(maxhead) & set(summary_header))
    return included_header


def checkrs(summary_data):
    from pkg_resources import resource_filename
    reference_file = resource_filename('cimr', 'data/annotation/variant_grch37_annotation_test.txt.gz')
    logging.info(f'using {reference_file} to check variant information.')
    reference = pandas.read_csv(reference_file, sep='\t', header=0, dtype={'chr':'str'})
    reference.columns = [x+'_reference' for x in reference.columns]
    rsnum_with_reference = summary_data.loc[summary_data['rsnum'].isin(reference['rs_id_dbSNP147_GRCh37p13_reference']),:]
    samples = rsnum_with_reference.sample(frac=0.2, replace=False)
    merged = samples.merge(reference, left_on='rsnum', right_on='rs_id_dbSNP147_GRCh37p13_reference', left_index=False, right_index=False, how='left')
    variant_nomatch = merged.loc[~(merged['variant_id']==merged['variant_id_reference'])]
    samplecount = len(merged.index)
    rsrefcount = len(rsnum_with_reference.index)
    nomatchcount = len(variant_nomatch.index)
    logging.info(f'out of {samplecount} sampled variants from {rsrefcount} total variants with RS IDs, {nomatchcount} variants do not match the reference.')
    return summary_data


def check_numeric(summary_data, col):
    """check for numeric columns"""
    from pandas.api.types import is_numeric_dtype
    try:
        if is_numeric_dtype(summary_data[col]):
            print(col+' is numeric')
            return summary_data
        else:
            numdata = (summary_data
                        .drop([col], axis=1)
                        .join(summary_data[col].apply(pandas.to_numeric, errors='coerce')))
            numcol = numdata[col].isnull().values().sum()
            logging.error(f'%s rows in %s are non-numeric'%(numcol,col,))
            return numdata
    except:
        logging.error(f'the format of %s is not testable.'%(col,))
        print(summary_data.head(n=2))
        sys.exit()


def readfile(filename):
    """read the input file as a pandas dataframe. check if empty"""
    filename = findfile(filename)
    summary_data = pandas.read_csv(filename, sep='\t', header=0)

    # check if empty and check header
    if not summary_data.empty:
        betaeffect = {'beta':'effect_size', 'se':'standard_error', 'pval':'pvalue'}
        summary_data.rename(columns=betaeffect, inplace=True)
        summary_header = summary_data.columns
        included_header = checkhead(summary_header)
    else:
        logging.error(f'no content in uploaded file {filename}.')    
        sys.exit()

    # check each column
    if 'variant_id' in included_header:
        summary_data = checkchrom(getpos(summary_data))
        logging.info(f'chromosome information is checked.')
    else:
        logging.error(f'variant_id column is not provided')
        pass
    if 'rsnum' in included_header:
        summary_data = checkrs(summary_data)
    else:
        logging.error(f'rsnum column is not provided.')
        pass
    if ('effect_size' and 'standard_error') in included_header:
        summary_data = check_numeric(check_numeric(summary_data, 'effect_size'), 'standard_error')
    else:
        pass
    return summary_data





