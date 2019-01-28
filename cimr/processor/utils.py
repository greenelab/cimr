#!/usr/bin/env python3

"""Utilities and common file checks used across different 
processor classes
(c) YoSon Park
"""

import sys
import pandas
import pathlib
import logging
import subprocess


class Infiler:
    """This is the cimr processor base class. functions regarding 
    automated checks for contributed summary statistics files are
    included in this class.

    Parameters
    ----------

    datatype: {'gwas', 'eqtl', 'sqtl', 'pqtl'}
    filename: name of the file to read in summary statistics
    genome_build: human genome reference id {'b37', 'b38'}


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
    GENOMEBUILDS = ('b37', 'b38')
    HEADERS = ['gene_id', 'rsnum', 'variant_id', 'pvalue', 
               'effect_size', 'standard_error', 'zscore', 'tss_distance', 
               'ma_samples', 'maf', 'chrom', 'pos', 'ref', 'alt', 
               'build', 'inc_allele', 'inc_afrq', 'imputation_status', 'n_cases', 'frequency'
               ]


    def __init__(self, datatype, filename, genome_build):
        if datatype not in self.DATATYPES:
            raise ValueError(' %s is not a valid datatype supported by cimr.' % datatype)
        if genome_build not in self.GENOMEBUILDS:
            raise ValueError(' %s is not a valid genome_build supported by cimr.' % genome_build)
        self.datatype = datatype
        self.filename = filename
        self.genome_build = genome_build
    

    def get_pos(self):
        """Check variant_id column and make 
        chrom, pos, ref, alt, build columns"""
        sumdata = self.summary_data
        temp = sumdata['variant_id'].str.split('_', expand=True)
        if not temp.empty:
            sumdata['chrom'] = temp[0]
            sumdata['pos'] = temp[1]
            sumdata['ref'] = temp[2]
            sumdata['alt'] = temp[3]
            sumdata['build'] = temp[4]
    

    def check_chrom(self, maxchrom=23):
        """Assumes chr+number
        - check for autosomal chromosomes
        - change if different from the specified format
        - discard non-autosomal chromosomes from main input
        """
        sumdata = self.summary_data
        chromdict = {str(i):'chr' + str(i) for i in range(1, maxchrom)}
        chromstr = ['chr' + str(i) for i in range(1, maxchrom)]
        chromint = [i for i in range(1, maxchrom)]
        chroms = sumdata['chrom'].drop_duplicates().values
        
        if len(chroms) > (maxchrom - 2) and len(chroms) < (maxchrom + 2):
            logging.info(f' there are {len(chroms)} chromosomes in the file provided.')
        elif len(chroms) <= (maxchrom-2):
            logging.warning(f' input file does not include {maxchrom} chromosome(s).')
            logging.warning(f' chromosome(s) included in the input file: %s'%(chroms,))
        else:
            logging.warning(f' input file more than {maxchrom - 1} chromosomes.')
            logging.warning(f' chromosome(s) included in the input file: %s'%(chroms,))

        if len(set(chroms) & set(chromstr)) > (maxchrom - 2):
            pass
        elif len(set(chroms) & set(chromint)) > (maxchrom - 2):
            sumdata['chrom'] = sumdata['chrom'].map(chromdict)
            sumdata = sumdata[sumdata['chrom'].isin(chromstr)]
            logging.info(f' chromosome ids have been updated.')
        else:
            logging.error(f' chromosome id needs to be checked.')

        remainder = list(set(chroms) - set(chromstr) - set(chromint))
        if len(remainder) > 0:
            logging.warning(f' chromosome(s) not used for analysis: %s'%(remainder,))


    def check_ref(self):
        from pkg_resources import resource_filename
        
        variant_reference_file = 'data/annotation/' + self.variant_reference_file
        reference_id = self.variant_reference_id
        reference_file = resource_filename(
            'cimr', variant_reference_file
            )
        logging.info(f' using {reference_file} to check variant information.')
        logging.info(f' rs id reference is {reference_id}')
        reference = pandas.read_csv(
            reference_file, sep='\t', header=0, dtype={'chr':'str'}
            )
        reference.columns = [x + '_reference' for x in reference.columns]
        sumdata = self.summary_data
        rsnum_with_reference = sumdata.loc[sumdata['rsnum'].isin(reference[reference_id + '_reference']),:]
        if not rsnum_with_reference.empty:
            samples = rsnum_with_reference.sample(frac=0.1, replace=False)
            merged = samples.merge(
                reference, left_on='rsnum', right_on=reference_id+'_reference', 
                left_index=False, right_index=False, how='left'
                )
            merged.drop_duplicates(inplace=True)
            variant_nomatch = merged.loc[~(merged['variant_id']==merged['variant_id_reference'])]
            samplecount = len(samples.index)
            rsrefcount = len(rsnum_with_reference.index)
            nomatchcount = len(variant_nomatch.index)
            logging.info(f' {samplecount} sampled variants from {rsrefcount} total variants with rs ids, {nomatchcount} variants do not match the reference.')
        else:
            logging.error(f' there are no matching rs ids.')
            pass


    def check_numeric(self, col):
        """Check for numeric columns"""
        from pandas.api.types import is_numeric_dtype
        try:
            if is_numeric_dtype(self.summary_data[col]):
                logging.info(f' {col} is numeric.')
            else:
                numdata = (self.summary_data
                            .drop([col], axis=1)
                            .join(self.summary_data[col].apply(pandas.to_numeric, errors='coerce')))
                numcol = numdata[col].isnull().values().sum()
                logging.error(f' %s rows in %s are non-numeric' % (numcol, col,))
                return numdata
        except:
            logging.error(f' the format of %s is not testable.' % (col,))
            print(self.summary_data.head(n=2))
            sys.exit()
    

    def check_probability(self, col):
        """Check whether probability is between 0 and 1"""
        if self.summary_data[col].between(0, 1, inclusive=True).any():
            logging.info(f' {str(col)} only contains values between 0 and 1.')
        else:
            logging.error(f' {str(col)} should only contain values between 0 and 1.')
        return 0
            

    def find_reference(self):
        """Find variant and gene references for map checking"""
        if self.genome_build == 'b37':
            self.gene_reference_file = 'gene_grch37_gencode_v29.txt.gz'
            self.variant_reference_file = 'variant_grch37_subset.txt.gz'
            self.variant_reference_id = 'rs_id_dbSNP147_GRCh37p13'
        else:
            self.gene_reference_file = 'gene_grch38_gencode_v26.txt.gz'
            self.variant_reference_file = 'variant_grch38_subset.txt.gz'
            self.variant_reference_id = 'rs_id_dbSNP150_GRCh38p7'
        return 0


    def read_file(self):
        """Read the input file as a pandas dataframe. check if empty"""
        if (pathlib.Path(self.filename).resolve()):
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
            sumdata = sumdata[self.included_header]
        else:
            logging.error(f' no content in uploaded file {self.filename}.')    
            sys.exit()
        
        self.find_reference()

        # check each column
        if 'variant_id' in self.included_header:
            self.get_pos()
            self.check_chrom()
            logging.info(f' chromosome information is checked.')
        else:
            logging.error(f' variant_id column is not provided')
            pass

        if 'rsnum' in self.included_header:
            self.check_ref()
        else:
            logging.error(f' rsnum column is not provided.')
            pass

        if 'effect_size' in self.included_header:
            self.check_numeric('effect_size') 
        else:
            logging.error(f' effect_size column is not provided.')
            pass

        if 'standard_error' in self.included_header:
            self.check_numeric('standard_error') 
        else:
            logging.error(f' standard_error column is not provided.')
            pass

        if 'pvalue' in self.included_header:
            self.check_numeric('pvalue')
            self.check_probability('pvalue')
        else:
            logging.error(f' pvalue column is not provided.')
            pass
        

    def write_file(self, outfile):
        """Write a checked file into a format used for cimr gene subprocess"""
        self.outfile = outfile
        
        try:
            self.summary_data.to_csv(self.outfile, header=True, index=False, sep='\t', na_reps='NA')
        except:
            logging.error(f' file {self.outfile} cannot be written.')
        return 0



class Integrator:
    """cimr integrator class connecting contributed data to cimr-adb

    Parameters:
    -----------

    filename: name of the file containing the data
    datatype = {'gwas', 'eqtl', 'tad'}
    can_be_public: boolean variable indicating whether the contributed data
                   can be release as a part of the public archive.
                   for cimr >= 0.2.3, can_be_public parameter will determine
                   the destination repository of the contributed data

    Notes:
    ------

    """


    def __init__(self, datatype, filename, can_be_public, genome_build):
        """File will be saved in cimr-adb"""
        self.datatype = datatype
        self.filename = filename
        self.can_be_public = can_be_public
        self.genome_build = genome_build


    def make_local_db(self, tempdir):
        """Temporarily download a local copy of the cimr-adb"""
        # TODO: change into public repo after testing
        self.tempdir = tempdir 
        try:
            clonercmd = 'git clone git@github.com:ypar/cimr-adb.git '+self.tempdir
            cloner = subprocess.Popen(
                clonercmd.split(), 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                shell=True, 
                executable='/bin/bash'
            )
            cloned = cloner.communicate()
            logging.info(f' {cloned}')
            logging.info(f' downloaded the existing database for Integrator')
        except OSError as e:
            print(e.errno)
            print(e.strerror)
        except:
            print(sys.exc_info()[0])
        return 0
    




