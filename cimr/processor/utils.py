#!/usr/bin/env python3

"""utilities and common file checks used across different 
processor classes
(c) YoSon Park
"""

import sys
import pandas
import pathlib
import logging


class Infiler:
    """cimr processor base class

    Parameters
    ----------

    datatype: {'gwas', 'eqtl', 'sqtl', 'pqtl'}
    filename: name of the file to read in summary statistics


    Notes:
    -------

    for files to be used for downstream functionalities,

    following are required:
    - gene_id : gene id. expected for eqtls
    - rsnum : rs id of the variant
    - variant_id : variant id in the following format
        chromosome_position_referenceallele_alternateallele_genomebuild
    - pval / pvalue: p-value of the beta estimate

    following are recommended:
    - beta / effect_size: beta coefficient estimate for the effect 
    of the variant on the gene 
    - se / standard_error: standard error of the beta
    - zscore: z-score if no beta/se information is present 
    (e.g. imputed summary statistic)

    following are optional:
    - tss_distance : distance to the transcription start site of the 
    gene_id
    - ma_samples : samples with minor alleles
    - maf : minor allele frequency
    - inc_allele : allele with respect to which variant's effect 
    sizes are estimated
    - inc_afrq : allele frequency of inc_allele

    following are given by parsing variant_id column:
    - chrom : chromosome id
    - pos: genomic position
    - ref : reference allele
    - alt : alternate allele
    - build : genomic build version number for the variant_id

    """

    DATATYPES = ('gwas', 'eqtl', 'sqtl', 'pqtl')
    HEADERS = ['gene_id', 'rsnum', 'variant_id', 'pvalue', 
               'effect_size', 'standard_error', 'zscore', 'tss_distance', 
               'ma_samples', 'maf', 'chrom', 'pos', 'ref', 'alt', 
               'build', 'inc_allele', 'inc_afrq'
               ]

    def __init__(self, datatype, filename):
        if datatype not in self.DATATYPES:
            raise ValueError(' %s is not a valid datatype supported by cimr.' % datatype)
        self.datatype = datatype
        self.filename = filename
    

    def getpos(self):
        """check variant_id column and make 
        chrom, pos, ref, alt, build columns"""
        sumdata = self.summary_data
        temp = sumdata['variant_id'].str.split('_', expand=True)
        if not temp.empty:
            sumdata['chrom'] = temp[0]
            sumdata['pos'] = temp[1]
            sumdata['ref'] = temp[2]
            sumdata['alt'] = temp[3]
            sumdata['build'] = temp[4]
    

    def checkchrom(self, maxchrom=23):
        """assumes chr+number
        - check for autosomal chromosomes
        - change if different from the specified format
        - discard non-autosomal chromosomes from main input
        """
        sumdata = self.summary_data
        chromdict = {str(i):'chr'+str(i) for i in range(1,maxchrom)}
        chromstr = ['chr'+str(i) for i in range(1,maxchrom)]
        chromint = [i for i in range(1,maxchrom)]
        chroms = sumdata['chrom'].drop_duplicates().values
        
        if len(chroms) > (maxchrom-2) and len(chroms) < (maxchrom+2):
            logging.info(f' there are {len(chroms)} chromosomes in the file provided.')
        elif len(chroms) <= (maxchrom-2):
            logging.warning(f' input file does not include {maxchrom} chromosome(s).')
            logging.warning(f' chromosome(s) included in the input file: %s'%(chroms,))
        else:
            logging.warning(f' input file more than {maxchrom-1} chromosomes.')
            logging.warning(f' chromosome(s) included in the input file: %s'%(chroms,))

        if len(set(chroms) & set(chromstr)) > (maxchrom-2):
            pass
        elif len(set(chroms) & set(chromint)) > (maxchrom-2):
            sumdata['chrom'] = sumdata['chrom'].map(chromdict)
            sumdata = sumdata[sumdata['chrom'].isin(chromstr)]
            logging.info(f' chromosome ids have been updated.')
        else:
            logging.error(f' chromosome id needs to be checked.')

        remainder = list(set(chroms) - set(chromstr) - set(chromint))
        if len(remainder) > 0:
            logging.warning(f' chromosome(s) not used for analysis: %s'%(remainder,))


    def checkrs(self):
        from pkg_resources import resource_filename
        reference_file = resource_filename(
            'cimr', 'data/annotation/variant_grch37_annotation_test.txt.gz'
            )
        logging.info(f'using {reference_file} to check variant information.')
        reference = pandas.read_csv(
            reference_file, sep='\t', header=0, dtype={'chr':'str'}
            )
        reference.columns = [x+'_reference' for x in reference.columns]
        sumdata = self.summary_data
        rsnum_with_reference = sumdata.loc[summary_data['rsnum'].isin(reference['rs_id_dbSNP147_GRCh37p13_reference']),:]
        samples = rsnum_with_reference.sample(frac=0.2, replace=False)
        merged = samples.merge(
            reference, left_on='rsnum', right_on='rs_id_dbSNP147_GRCh37p13_reference', 
            left_index=False, right_index=False, how='left'
            )
        variant_nomatch = merged.loc[~(merged['variant_id']==merged['variant_id_reference'])]
        samplecount = len(merged.index)
        rsrefcount = len(rsnum_with_reference.index)
        nomatchcount = len(variant_nomatch.index)
        logging.info(f' {samplecount} sampled variants from {rsrefcount} total variants with RS IDs, {nomatchcount} variants do not match the reference.')


    def check_numeric(self, col):
        """check for numeric columns"""
        from pandas.api.types import is_numeric_dtype
        try:
            if is_numeric_dtype(self.summary_data[col]):
                print(col+' is numeric')
                return self
            else:
                numdata = (self.summary_data
                            .drop([col], axis=1)
                            .join(self.summary_data[col].apply(pandas.to_numeric, errors='coerce')))
                numcol = numdata[col].isnull().values().sum()
                logging.error(f' %s rows in %s are non-numeric'%(numcol,col,))
                return numdata
        except:
            logging.error(f' the format of %s is not testable.'%(col,))
            print(self.summary_data.head(n=2))
            sys.exit()


    def readfile(self):
        """read the input file as a pandas dataframe. check if empty"""
        if(pathlib.Path(self.filename).resolve()):
            self.filename = str(self.filename)
            logging.info(f' processing {self.filename}.')
        else:
            logging.error(f' no file {self.filename} found for processing.')
            sys.exit()

        self.summary_data = pandas.read_csv(self.filename, sep='\t', header=0)
        sumdata = self.summary_data

        # check if empty and check header
        if not sumdata.empty:
            betaeffect = {'beta':'effect_size', 'se':'standard_error', 'pval':'pvalue'}
            sumdata.rename(columns=betaeffect, inplace=True)
            self.included_header = list(set(self.HEADERS) & set(self.summary_data.columns))
        else:
            logging.error(f'no content in uploaded file {self.filename}.')    
            sys.exit()

        # check each column
        if 'variant_id' in self.included_header:
            self.getpos()
            self.checkchrom()
            logging.info(f' chromosome information is checked.')
        else:
            logging.error(f' variant_id column is not provided')
            pass

        if 'rsnum' in self.included_header:
            self.checkrs()
        else:
            logging.error(f' rsnum column is not provided.')
            pass

        if ('effect_size' and 'standard_error') in self.included_header:
            self.check_numeric('effect_size') 
            self.check_numeric('standard_error')
        else:
            pass






