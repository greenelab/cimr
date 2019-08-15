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

from pandas.api.types import is_numeric_dtype

from .query import Querier

from ..defaults import DATA_TYPES
from ..defaults import GENOME_BUILDS
from ..defaults import HEADER
from ..defaults import MAXCHROM


def set_chrom_dict():
    """Make a dictionary to standardize chromosome IDs in input files."""
    chrom_dict = {
        str(i):'chr' + str(i) for i in range(1, MAXCHROM)
    }
    chrom_dict.update({
        'X':'chr23', 
        'Y':'chr24', 
        'M':'chr25', 
        'MT':'chr25', 
        'chrX':'chr23', 
        'chrY':'chr24', 
        'chrM':'chr25', 
        'chrMT':'chr25'
    })
    return chrom_dict, MAXCHROM


def find_file(file_name):
    """Check if a file exists and exit if not."""
    if (pathlib.Path(file_name).resolve()):
        file_name = str(file_name)
        logging.info(f' processing {file_name}.')
        return file_name
    else:
        logging.error(f' no file {file_name} found for processing.')
        sys.exit()


def check_numeric(data, col):
    """Check for numeric columns"""
    from pandas.api.types import is_numeric_dtype
    try:
        if is_numeric_dtype(data[col]):
            logging.info(f' {col} is numeric.')
        else:
            numdata = (data
                        .drop([col], axis=1)
                        .join(data[col].apply(pandas.to_numeric, errors='coerce'))
                        )
            numcol = numdata[col].isnull().values().sum()
            logging.error(f' %s rows in %s are non-numeric' % (numcol, col,))
            return numdata
    except:
        logging.error(f' the format of %s is not testable.' % (col,))
        print(data.head(n=2))
        sys.exit(1)


def make_int(col):
    """Make sure the column values are represented as int"""
    return col.astype(int)


class Infiler:
    """This is the cimr processor base class. functions regarding 
    automated checks for contributed summary statistics files are
    included in this class.

    Parameters
    ----------

    data_type: {'gwas', 'eqtl', 'sqtl', 'pqtl', 'twas', 'tad'}
    file_name: name of the file to read in summary statistics
    genome_build: human genome reference id {'b37', 'b38'}


    Notes:
    -------

    for files to be used for downstream functionalities,

    following are required :
    - feature_id : e.g. gene id. expected for eqtls
    - rsnum : rs id of the variant
    - variant_id : variant id in the following format
        chromosome_position_referenceallele_alternateallele_genomebuild
    - pval / pvalue: p-value of the beta estimate
    - beta / effect_size : beta coefficient estimate for the effect 
        of the variant on the gene 
    - se / standard_error : standard error of the beta
    - zscore: z-score if no beta/se information is present 
        (e.g. imputed summary statistic)

    following are optional :
    - tss_distance : distance to the transcription start site of the 
        gene_id
    - ma_samples : samples with minor alleles
    - maf : minor allele frequency
    - effect_allele : allele with respect to which variant's effect 
        sizes are estimated
    - non_effect_allele : non-effect-allele
    - imputation_status : imputation status of the variant. 'original' to
        indicate that the original GWAS contained the variant. 'imputed' to
        indicate that the data submitter has done (additional) imputation 
        to the initial release of the data.
    - sample_size : sample size for the variant underlying the effect_size 
        and pvalue estimates. If not available per-variant, the total sample
        size used for the study is accepted.
    - n_cases : for binary phenotypic traits, n_cases are necessary for 
        downstream analyses such as coloc colocalization.
    - frequency : allele frequency of the effect_allele

    following are given by parsing variant_id column:
    - chrom : chromosome id
    - pos: genomic position
    - ref : reference allele
    - alt : alternate allele
    - build : genomic build version number for the variant_id

    """

    def __init__(self, 
                 data_type, 
                 file_name, 
                 genome_build, 
                 update_rsid, 
                 outfile, 
                 chunksize,
                 columnset):

        if data_type not in DATA_TYPES:
            raise ValueError(' %s is not a valid data_type supported' % data_type)
        if genome_build not in GENOME_BUILDS:
            raise ValueError(' %s is not a valid genome_build supported' % genome_build)
        self.data_type = data_type
        self.file_name = file_name
        self.genome_build = genome_build
        self.update_rsid = update_rsid
        self.outfile = outfile
        self.chunksize = chunksize
        self.columnset = columnset
    

    def get_pos(self):
        """Check variant_id column and make 
        chrom, pos, ref, alt, build columns"""
        # check first element to find delimiter in variant_id
        variant_ids = self.summary_data['variant_id']
        if '_' in variant_ids[0]:
            pass
        elif '-' in variant_ids[0]:
            variant_ids = variant_ids.str.replace('-', '_')
            logging.info(f' variant_id delimiter changed from \'-\' to \'_\'')
        elif ':' in variant_ids[0]:
            variant_ids = variant_ids.str.replace(':', '_')
            logging.info(f' variant_id delimiter changed from \':\' to \'_\'')
        else:
            logging.error(f' unknown delimiter used in variant_id')
            sys.exit(1)

        temp = variant_ids.str.split('_', expand=True)

        if not temp.empty:
            self.summary_data['chrom'] = temp[0]
            self.summary_data['pos'] = temp[1]
            self.summary_data['ref'] = temp[2]
            self.summary_data['alt'] = temp[3]
            if len(temp.columns) == 5:
                self.summary_data['build'] = temp[4]
            else:
                logging.info(f' updating variant_id to include build')
                self.summary_data['build'] = self.genome_build
                self.summary_data['variant_id'] = variant_ids + '_' + self.genome_build
    

    def make_variant_id(self):
        """(Re)make variant_id with updated chrom ID, build #, etc."""
        self.summary_data['variant_id'] = self.summary_data['chrom'] \
            + '_' \
            + self.summary_data['pos'] \
            + '_' \
            + self.summary_data['ref'] \
            + '_' \
            + self.summary_data['alt'] \
            + '_' \
            + self.summary_data['build']


    def check_chrom(self):
        """Assumes chr+number
        * check for autosomal chromosomes
        * change if different from the specified format
        * discard non-autosomal chromosomes from main input
        """
        sumdata = self.summary_data
        chrom_dict, maxchrom = set_chrom_dict()
        chrom_str = list(chrom_dict.values())[0:maxchrom]
        chrom_int = list(chrom_dict.keys())[0:maxchrom]
        chroms = sumdata['chrom'].drop_duplicates().values
        
        if len(chroms) > (maxchrom - 2) and len(chroms) < (maxchrom + 2):
            logging.info(f' there are {len(chroms)} chromosomes.')
        elif len(chroms) <= (maxchrom - 2):
            logging.warning(f' chromosome(s) included: %s'%(chroms,))
        else:
            logging.warning(f' input file more than {maxchrom - 1} chromosomes.')
            logging.warning(f' chromosome(s) included: %s'%(chroms,))

        if len(set(chroms) & set(chrom_str)) > (maxchrom - 2):
            pass
        elif len(set(chroms) & set(chrom_int)) > (maxchrom - 2):
            sumdata['chrom'] = sumdata['chrom'].map(
                chrom_dict, na_action='ignore'
            ).fillna(sumdata['chrom'])
            sumdata = sumdata[sumdata['chrom'].isin(chrom_str)]
            logging.info(f' chromosome ids have been updated.')
        else:
            logging.error(f' chromosome id needs to be checked.')

        remainder = list(set(chroms) - set(chrom_str) - set(chrom_int))
        if len(remainder) > 0:
            logging.warning(f' chromosome(s) not used: %s'%(remainder,))


    def check_ref(self):
        from pkg_resources import resource_filename
        
        var_ref = 'data/annotation/' + self.var_ref
        ref_id = self.var_ref_id
        ref_file = resource_filename(
            'cimr', var_ref
            )
        logging.info(f' using {ref_file} to check variants.')
        logging.info(f' rs id reference is {ref_id}')
        refdf = pandas.read_csv(
            ref_file, sep='\t', header=0, dtype={'chr':'str'}
            )
        refdf.columns = [x + '_reference' for x in refdf.columns]
        sumdata = self.summary_data
        checked_ref = sumdata['rsnum'].isin(refdf[ref_id + '_reference'])
        rsnum_ref = sumdata.loc[checked_ref,:]
        if not rsnum_ref.empty:
            samples = rsnum_ref.sample(frac=0.1, replace=False)
            merged = samples.merge(
                refdf, left_on='rsnum', right_on=ref_id+'_reference', 
                left_index=False, right_index=False, how='left'
            )
            merged.drop_duplicates(inplace=True)
            variant_nomatch = merged.loc[~(merged['variant_id']==merged['variant_id_reference'])]
            samplecount = len(samples.index)
            rsrefcount = len(rsnum_ref.index)
            nomatchcount = len(variant_nomatch.index)
            logging.info(f' {samplecount} sampled from {rsrefcount} total variants with rs ids,')
            logging.info(f' {nomatchcount} variants do not match the reference.')
        else:
            logging.error(f' there are no matching rs ids.')
            pass
    

    def check_probability(self, col):
        """Check whether probability is between 0 and 1"""
        if self.summary_data[col].between(0, 1, inclusive=True).any():
            logging.info(f' {str(col)} only contains values between 0 and 1.')
        else:
            logging.error(f' {str(col)} should only contain values between 0 and 1.')
            

    def find_reference(self):
        """Find variant and gene references for map checking"""
        if self.genome_build == 'b37':
            self.gene_ref = 'gene_grch37_gencode_v26.tsv.gz'
            self.var_ref = 'variant_grch37_subset.txt.gz'
            self.var_ref_id = 'rs_id_dbSNP147_GRCh37p13'
        else:
            self.gene_ref = 'gene_grch38_gencode_v29.tsv.gz'
            self.var_ref = 'variant_grch38_subset.txt.gz'
            self.var_ref_id = 'rs_id_dbSNP150_GRCh38p7'
    

    def trim_ensembl(self):
        """Make the ensembl gene id list to be agnostic to 
        subversions for queries
        
        e.g. ENSG00000143614.7 -> ENSG00000143614
        """
        ensemblid = self.summary_data['feature_id'].str.split('.').str[0]
        self.summary_data['feature_id'] = ensemblid


    def list_features(self):
        """Find the list of features (e.g. genes)."""
        logging.debug(f' {self.included_header}')
        logging.debug(f' {self.summary_data.head(2)}')
        if 'feature_id' in self.included_header:
            if self.summary_data['feature_id'][0].startswith('ENSG'):
                self.trim_ensembl()

        if 'ensemblgene' in self.summary_data.columns:
            return self.summary_data.ensemblgene
        elif 'feature_id' in self.summary_data.columns:
            return self.summary_data.feature_id
        else:
            logging.error(f' feature_id column is not provided.')
            sys.exit(1)
    

    def append_gene_cols(self, gene_df):
        """Add gene annotation columns to dataframe"""
        self.summary_data = pandas.concat([self.summary_data, gene_df], axis=1)
    

    def fill_effect_allele(self):
        """Fill NA in effect_allele column with inc_allele if 
        inc_allele present
        """
        self.summary_data['effect_allele'].fillna(self.summary_data['inc_allele'])
        

    def make_int(self, colname):
        """Make int columns int"""
        self.summary_data[colname] = self.summary_data[colname].astype(int)


    def check_file(self, summary_data):
        """Check different columns for dtype, remove missing rows and 
        standardize format to be used for analyses
        """
        
        self.find_reference()
        summary_data.reset_index(inplace=True, drop=True)
        self.summary_data = summary_data

        if 'variant_id' in self.included_header:
            self.get_pos()
            self.check_chrom()
            logging.info(f' chromosome information is checked.')
            self.make_variant_id()
            logging.info(f' variant_id has been standardized.')
        else:
            logging.error(f' variant_id column is not provided')
            sys.exit(1)

        if 'rsnum' in self.included_header:
            if self.update_rsid:
                self.check_ref()
        else:
            logging.warning(f' rsnum column is not provided.')
            
        if 'ma_samples' in self.included_header:
            self.make_int('ma_samples')
        
        if 'ma_count' in self.included_header:
            self.make_int('ma_count')
        
        if 'sample_size' in self.included_header:
            self.make_int('sample_size')

        if 'inc_allele' in self.included_header:
            self.fill_effect_allele()
        else:
            logging.debug(f' inc_allele column is not available')
        
        columns_to_drop = [
            'chrom', 'pos', 'ref', 'alt', 'chromosome', 'build',
            'position', 'inc_allele', 'variant_chrom', 'variant_pos'
        ]

        for colname in columns_to_drop:
            if colname in self.summary_data.columns:
                self.summary_data.drop(
                    colname, 
                    axis=1, 
                    inplace=True
                )

        if 'effect_size' in self.included_header:
            check_numeric(self.summary_data, 'effect_size') 
        else:
            logging.error(f' effect_size column is not provided.')
            sys.exit(1)

        if 'standard_error' in self.included_header:
            check_numeric(self.summary_data, 'standard_error') 
        else:
            logging.error(f' standard_error column is not provided.')
            sys.exit(1)

        if 'pvalue' in self.included_header:
            check_numeric(self.summary_data, 'pvalue')
            self.check_probability('pvalue')
        else:
            logging.error(f' pvalue column is not provided.')
            sys.exit(1)

        if 'pvalue_perm' in self.included_header:
            check_numeric(self.summary_data, 'pvalue_perm')
            self.check_probability('pvalue_perm')
        else:
            logging.debug(f' pvalue_perm column is not provided.')

        if 'fdr' in self.included_header:
            check_numeric(self.summary_data, 'fdr')
            self.check_probability('fdr')
        else:
            logging.debug(f' fdr column is not provided.')

        if 'qvalue' in self.included_header:
            check_numeric(self.summary_data, 'qvalue')
        else:
            logging.debug(f' qvalue column is not provided.')


    def write_header(self):
        """Write data to file with header"""
        if isinstance(self.summary_data, pandas.DataFrame):
            self.summary_data.to_csv(
                str(self.outfile), 
                header=True, 
                index=False, 
                sep='\t', 
                na_rep='NA',
                compression='gzip',
                float_format='%.5f',
                mode='w'
            )


    def write_file(self):
        """Write data to file without header"""
        if isinstance(self.summary_data, pandas.DataFrame):
            self.summary_data.to_csv(
                str(self.outfile), 
                header=False, 
                index=False, 
                sep='\t', 
                na_rep='NA',
                compression='gzip',
                float_format='%.5f',
                mode='a'
            )
    
    def check_h5_outfile(self):
        """Check file extension for h5 output"""
        if not str(self.outfile).endswith('.h5.bz2'):
            self.outfile = str(self.outfile) + '.h5.bz2'


    def write_h5_header(self):
        """Write data as PyTable in an h5 file system with header"""
        if isinstance(self.summary_data, pandas.DataFrame):
            self.summary_data.to_hdf(
                str(self.outfile), 
                header=True, 
                index=False, 
                format='table', 
                complevel=9, 
                complib='bzip2', 
                key=self.data_type, 
                mode='w'
            )

    def write_h5_file(self):
        """Append data as PyTable in an h5 file system"""
        if isinstance(self.summary_data, pandas.DataFrame):
            self.summary_data.to_hdf(
                str(self.outfile), 
                header=True, 
                index=False, 
                format='table', 
                complevel=9, 
                complib='bzip2', 
                key=self.data_type, 
                mode='a'
            )


    def rename_columns(self, dataframe):
        """Given a pandas dataframe and a dictionary, rename columns
        Primarily used for yaml column-set-based renaming in cimr.
        """
        logging.info(f' renaming columns based on provided info')
        dataframe.rename(self.columnset, axis=1, inplace=True)
    

    def call_querier(self, genes):
        queried = Querier(genes)
        queried.form_query()
        self.append_gene_cols(queried.list_queried())
    

    def map_features(self, features):
        """Find the type of feature_id, map to gene symbols, 
        if not available.
        """
        if features.str.startswith('ENSG').all(skipna=True):
            self.feature_type = 'ensemblgene'
        elif features.str.startswith('ENST').all(skipna=True):
            self.feature_type = 'ensembltranscript'
        elif features.str.startswith('ENSP').all(skipna=True):
            self.feature_type = 'ensemblprotein'
        elif is_numeric_dtype(features):
            self.feature_type = 'entrezgene'
        else:
            logging.error(f' feature reference cannot be determined.')
            sys.exit(1)
        
        if not 'feature_name' in self.included_header:
            gene_annot = pandas.read_csv(
                'cimr/data/annotation/' + self.gene_ref, 
                sep='\t', 
                header=0
            )
            cols = ['feature_name', self.feature_type, 'feature_type']
            gene_annot = gene_annot[cols]
            gene_annot.rename(columns={self.feature_type:'feature_id'}, inplace=True)

            logging.debug(f' current dataframe: ')
            logging.debug(f' {self.summary_data.head(2)}')
            logging.debug(f' selected annotations: ')
            logging.debug(f' {gene_annot.head(2)}')
            
            self.summary_data = self.summary_data.merge(
                gene_annot, 
                on='feature_id', 
                how='left',
                left_index=False,
                right_index=False
            )
            logging.debug(f' dataframe has been annotated.')
            logging.debug(f' {self.summary_data.head(2)}')
    

    def read_file(self):
        """Read the input file as a pandas dataframe. check if empty"""
        self.file_name = find_file(self.file_name)

        chunks = pandas.read_csv(
            self.file_name, 
            sep='\t', 
            header=0, 
            iterator=True,
            chunksize=self.chunksize
        )

        chunkcount = 0

        for chunk in chunks:
            chunk.reset_index(drop=True, inplace=True)
            logging.info(f' processing input chunk {chunkcount}')
            # check if empty and check header
            if not chunk.empty:
              
                chunk.reset_index(drop=True, inplace=True)
                if self.columnset:
                    self.rename_columns(chunk)

                self.included_header = list(
                    set(HEADER) & set(chunk.columns)
                )
                self.check_file(chunk)
                logging.debug(f' processing data type {self.data_type}')

                if self.data_type == 'eqtl':
                    features = self.list_features()
                    # 413 error
                    # self.call_querier(features)
                    self.map_features(features)
                
                logging.info(f' dropping duplicate columns...')
                dropcols = self.summary_data.columns.duplicated()
                self.summary_data = self.summary_data.loc[:, ~dropcols]

                logging.info(f' writing processed data...')

                if chunkcount == 0:
                    self.write_header()
                elif chunkcount > 0:
                    self.write_file()
                else:
                    logging.error(f' file {self.outfile} cannot be written.')
                    sys.exit(1)
                
            else:
                logging.error(f' no content in {self.file_name}.')
                sys.exit(1)
            
            chunkcount += 1
        

class Integrator:
    """cimr integrator class connecting contributed data to cimr-d

    Parameters:
    -----------

    data_type: {'gwas', 'eqtl', 'sqtl', 'pqtl', 'twas', 'tad'}
    file_name: name of the file containing the data
    data_type = defaults.DATA_TYPES
    can_be_public: boolean variable indicating whether the contributed 
                   data can be released as a part of the public archive.
                   for cimr >= 0.2.3, can_be_public parameter will 
                   determine the destination repository

    Notes:
    ------
    Integrator is primarily used by creating a PR with new data in cimr-d.

    """


    def __init__(self, data_type, file_name, can_be_public, genome_build):
        """File will be saved in cimr-d"""
        self.data_type = data_type
        self.file_name = file_name
        self.can_be_public = can_be_public
        self.genome_build = genome_build


    def make_local_db(self, tempdir):
        """Temporarily download a local copy of the cimr-d"""
        # TODO: change into public repo after testing
        self.tempdir = tempdir 
        try:
            clonercmd = 'git clone git@github.com:greenelab/cimr-d.git ' + self.tempdir
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


