#!/usr/bin/env python3

"""cimr processor submodule main classes

(c) YoSon Park
"""

import sys
import pandas
import pathlib
import logging
import subprocess

from pandas.api.types import is_numeric_dtype

# Querier class for feature_id mapping
from .query import Querier

# utilities to convert values
from .convertibles import (get_effect_direction, get_z,
    convert_p_to_z, convert_z_to_p, convert_or_to_beta, estimate_se)

# utilities
from .utils import (set_chrom_dict, find_file, intersect_set,
    check_numeric, check_probability)

# default values
from ..defaults import (COMPRESSION_EXTENSION, ANNOTURL,
    DATA_TYPES, GENOME_BUILDS, VAR_COMPONENTS, MAXCHROM,
    HEADER, REQ_COLUMNS, NUMERIC_COLUMNS, PROB_COLUMNS, INT_COLUMNS)



class Infiler:
    """This is the cimr processor base class. functions regarding
    automated checks for contributed summary statistics files are
    included in this class.

    ------------

    data_type: {'gwas', 'eqtl', 'sqtl', 'pqtl', 'twas', 'tad'}
    file_name: name of the file to read in summary statistics
    genome_build: human genome reference id {'b37', 'b38'}
    data: pandas DataFrame of {'chunksize'}
    update_rsid: whether to update rsid, currently deprecated
    outfile: string variable for outfile name
    chunksize: chunksize to be read at a time for pandas.read_csv()
    columnset: a dictionary of columns to be renamed


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
        self.data = pandas.DataFrame()
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
        variant_ids = self.data['variant_id']
        if '_' in variant_ids[0]:
            pass
        if '-' in variant_ids[0]:
            variant_ids = variant_ids.str.replace('-', '_')
            logging.info(f' variant_id delimiter changed from \'-\' to \'_\'')
        if ':' in variant_ids[0]:
            variant_ids = variant_ids.str.replace(':', '_')
            logging.info(f' variant_id delimiter changed from \':\' to \'_\'')
        # else:
        #     logging.error(f' unknown delimiter used in variant_id')
        #     sys.exit(1)

        temp = variant_ids.str.split('_', expand=True)

        if not temp.empty:
            self.data['chrom'] = temp[0]
            self.data['pos'] = temp[1]
            self.data['ref'] = temp[2]
            self.data['alt'] = temp[3]
            if len(temp.columns) == 5:
                self.data['build'] = temp[4]
            else:
                logging.info(f' updating variant_id to include build')
                self.data['build'] = self.genome_build
                self.data['variant_id'] = variant_ids + '_' + self.genome_build


    def make_variant_id(self):
        """(Re)make variant_id with updated chrom ID, build #, etc."""
        logging.debug(f' standardizing allele codes...')
        self.data['ref'] = self.data['ref'].str.upper()
        self.data['alt'] = self.data['alt'].str.upper()

        logging.debug(f' making a new variant_id column.')

        self.data['variant_id'] = self.data['chrom'].astype(str) \
            + '_' \
            + self.data['pos'].astype(str) \
            + '_' \
            + self.data['ref'].astype(str) \
            + '_' \
            + self.data['alt'].astype(str) \
            + '_' \
            + self.data['build'].astype(str)
        logging.debug(f' variant_id column verified.')
        logging.info(f' variant_id has been standardized.')


    def check_chrom(self):
        """Assumes chr+number
        * change if different from the specified format
        * check if integer turned into float
          e.g. 1 -> 1.0
        """
        chrom_dict, maxchrom = set_chrom_dict()
        chrom_str = list(chrom_dict.values())[0:maxchrom]
        chrom_int = list(chrom_dict.keys())[0:maxchrom]
        # for files with integer value chromosomes
        # (i.e. 1 in stead of 'chr1'),
        # sometimes values are turned into float, which will not
        # work with default chrom_dict.
        # creating a chrom_flt to match these, check if int, then
        # standardize into str 'chr' + int format
        chrom_flt = [float(x) for x in range(1, maxchrom)]

        chroms = self.data['chrom'].drop_duplicates().values

        if len(set(chroms) & set(chrom_flt)) > 0:
            self.data['chrom'] = self.data['chrom'].astype(int)
            self.data['chrom'] = self.data['chrom'].astype(str)
            chroms = self.data['chrom'].drop_duplicates().values
            logging.debug(f' {chroms}')

        if len(chroms) > (maxchrom - 2) and len(chroms) < (maxchrom + 2):
            logging.info(f' there are {len(chroms)} chromosomes.')
        elif len(chroms) <= (maxchrom - 2):
            logging.info(f' chromosome(s) included: %s'%(chroms,))
        else:
            logging.warning(f' input file more than {maxchrom - 1} chromosomes.')
            logging.warning(f' chromosome(s) included: %s'%(chroms,))

        if len(set(chroms) & set(chrom_str)) > (maxchrom - 2):
            pass

        if len(set(chroms) & set(chrom_int)) > 0:
            self.data['chrom'] = self.data['chrom'].map(
                chrom_dict, na_action='ignore'
            ).fillna(self.data['chrom'])
            logging.info(f' chromosome ids have been updated.')
        else:
            logging.warning(f' chromosome id needs to be checked.')

        remainder = list(set(chroms) - set(chrom_str) - set(chrom_int))
        if len(remainder) > 0:
            logging.warning(f' chromosome(s) not used: %s'%(remainder,))

        logging.info(f' finished checking chromosome ids.')


    def check_ref(self):
        """Compare a random subset of variants to GRCh37 or GRCh38 references."""
        from pkg_resources import resource_filename

        var_ref = self.var_ref
        ref_id = self.var_ref_id
        ref_file = resource_filename(
            'cimr', var_ref
            )
        logging.info(f' using {ref_file} to check variants.')
        logging.info(f' rs id reference is {ref_id}.')
        refdf = pandas.read_csv(
            ref_file, sep='\t', header=0, dtype={'chr':'str'}
        )
        refdf.columns = [x + '_reference' for x in refdf.columns]
        sumdata = self.data.copy()
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


    def find_reference(self):
        """Find variant and gene references for map checking"""
        if self.genome_build == 'b37':
            self.gene_ref = ANNOTURL + 'gene_grch37_gencode_v26.tsv.gz'
            self.var_ref = ANNOTURL + 'variant_grch37_subset.txt.gz'
            self.var_ref_id = 'rs_id_dbSNP147_GRCh37p13'
        else:
            self.gene_ref = ANNOTURL + 'gene_grch38_gencode_v29.tsv.gz'
            self.var_ref = ANNOTURL + 'variant_grch38_subset.txt.gz'
            self.var_ref_id = 'rs_id_dbSNP150_GRCh38p7'


    def trim_ensembl(self):
        """Make the ensembl gene id list to be agnostic to
        subversions for queries

        e.g. ENSG00000143614.7 -> ENSG00000143614
        """
        ensemblid = self.data['feature_id'].str.split('.').str[0]
        self.data['feature_id'] = ensemblid
        logging.info(f' ensembl id has been truncated for database queries.')


    def list_features(self):
        """Find the list of features (e.g. genes)."""
        logging.debug(f' {self.included_header}')
        logging.debug(f' {self.data.head(2)}')
        if 'feature_id' in self.included_header:
            if self.data['feature_id'][0].startswith('ENS'):
                self.trim_ensembl()

        if 'ensemblgene' in self.data.columns:
            return self.data.ensemblgene
        elif 'feature_id' in self.data.columns:
            return self.data.feature_id
        else:
            logging.error(f' feature_id column is not provided.')
            sys.exit(1)


    def append_gene_cols(self, gene_df):
        """Add gene annotation columns to dataframe"""
        self.data = pandas.concat(
            [self.data, gene_df],
            axis=1
        )


    def fill_effect_allele(self):
        """Fill NA in effect_allele column with inc_allele if
        inc_allele present
        """
        self.data['effect_allele'].fillna(
            self.data['inc_allele']
        )


    def make_composite_id(self):
        """Files may be missing the standardized variant_id
        but may contain chrom, pos, ref and alt columns to make a
        new variant_id column. This option can be activated
        by indicating all variant_id component fields in the yaml.
        """
        logging.debug(f' processing information needed to make variant_id.')

        if set(VAR_COMPONENTS).issubset(set(self.included_header)):
            logging.debug(f' checking {VAR_COMPONENTS} for missing values.')
            self.data.dropna(
                subset=VAR_COMPONENTS, inplace=True
            )
            logging.debug(f' rows with missing {VAR_COMPONENTS} are dropped.')

            logging.debug(f' making variant_id using {VAR_COMPONENTS}.')
            self.data['chrom'] = self.data['variant_chrom']
            self.check_chrom()
            logging.info(f' verifying variant positions are int values.')
            self.data['variant_pos'] = self.data['variant_pos'].astype(int)
            logging.debug(f' changing verified values to str.')
            self.data['variant_pos'] = self.data['variant_pos'].astype(str)

            self.data['pos'] = self.data['variant_pos']
            self.data['ref'] = self.data['non_effect_allele']
            self.data['alt'] = self.data['effect_allele']
            self.data['build'] = self.genome_build

            # make variant_id from components
            self.make_variant_id()
            # recall included_header from columns including variant_id
            self.included_header = intersect_set(
                HEADER, self.data.columns
            )

        else:
            logging.error(f' missing columns necessary to make variant_id.')
            sys.exit(1)


    def check_data(self, data):
        """Check different columns for dtype, remove missing rows and
        standardize format to be used for analyses
        """
        logging.debug(f' HEADER: {HEADER}.')
        logging.debug(f' columns: {data.columns}.')
        self.included_header = intersect_set(
            HEADER, data.columns
        )
        logging.debug(f' included header overlapping cimr default set: {self.included_header}.')

        self.find_reference()
        data.reset_index(inplace=True, drop=True)
        self.data = data.copy()

        if 'variant_id' not in self.data.columns:
            self.make_composite_id()
            logging.debug(f' data.head(2): {self.data.head(2)}')
        elif 'variant_id' in self.included_header:
            self.get_pos()
            self.check_chrom()
            self.make_variant_id()
        else:
            logging.error(f' variant_id column is not provided.')
            sys.exit(1)

        if 'rsnum' in self.data.columns:
            if self.update_rsid:
                self.check_ref()
        else:
            logging.warning(f' rsnum column is not provided.')

        for col in INT_COLUMNS:
            if col in self.data.columns:
                self.data[col] = self.data[col].astype(int)

        if 'inc_allele' in self.data.columns:
            self.fill_effect_allele()
        else:
            logging.debug(f' inc_allele column is not available.')

        columns_to_drop = [
            'chrom', 'pos', 'ref', 'alt', 'chromosome', 'build',
            'position', 'inc_allele', 'variant_chrom', 'variant_pos'
        ]

        for col in columns_to_drop:
            if col in self.data.columns:
                self.data.drop(
                    col,
                    axis=1,
                    inplace=True
                )

        self.data = estimate_se(self.data)
        self.data = convert_z_to_p(self.data)
        self.data = convert_p_to_z(self.data)
        self.data = convert_or_to_beta(self.data)
        self.data = get_z(self.data)

        for col in REQ_COLUMNS:
            if col not in self.data.columns:
                logging.error(f' {col} is required.')
                sys.exit(1)

        for col in NUMERIC_COLUMNS:
            if col in self.data.columns:
                self.data = check_numeric(self.data, col)

        for col in PROB_COLUMNS:
            if col in self.data.columns:
                check_probability(self.data, col)


    def write_header(self):
        """Write data to file with header"""
        if isinstance(self.data, pandas.DataFrame):
            self.data.to_csv(
                str(self.outfile),
                header=True,
                index=False,
                sep='\t',
                na_rep='NA',
                compression='gzip',
                float_format='%.6f',
                mode='w'
            )


    def write_file(self):
        """Write data to file without header"""
        if isinstance(self.data, pandas.DataFrame):
            self.data.to_csv(
                str(self.outfile),
                header=False,
                index=False,
                sep='\t',
                na_rep='NA',
                compression='gzip',
                float_format='%.6f',
                mode='a'
            )

    def check_h5_outfile(self):
        """Check file extension for h5 output"""
        if not str(self.outfile).endswith('.h5.bz2'):
            self.outfile = str(self.outfile) + '.h5.bz2'


    def write_h5_header(self):
        """Write data as PyTable in an h5 file system with header"""
        if isinstance(self.data, pandas.DataFrame):
            self.data.to_hdf(
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
        if isinstance(self.data, pandas.DataFrame):
            self.data.to_hdf(
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
        logging.info(f' renaming columns based on column dictionary.')
        dataframe.rename(self.columnset, axis=1, inplace=True)
        logging.debug(f' renamed data.head(2): {dataframe.head(2)}')


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

        logging.debug(f' reading {self.gene_ref}.')

        if not 'feature_name' in self.included_header:
            gene_annot = pandas.read_csv(
                self.gene_ref,
                sep='\t',
                header=0
            )
            cols = ['feature_name', self.feature_type, 'feature_type']
            gene_annot = gene_annot[cols]
            gene_annot.rename(columns={self.feature_type:'feature_id'}, inplace=True)

            logging.debug(f' current dataframe: {self.data.head(2)}')
            logging.debug(f' selected annotations: {gene_annot.head(2)}')

            self.data = self.data.merge(
                gene_annot,
                on='feature_id',
                how='left',
                left_index=False,
                right_index=False
            )
            logging.debug(f' dropping duplicates, if any...')
            self.data.drop_duplicates(inplace=True)
            self.data.reset_index(inplace=True, drop=True)
            logging.debug(f' dataframe has been annotated.')
            logging.debug(f' {self.data.head(2)}')


    def order_columns(self):
        """Make sure the final output has the required columns
        listed first."""
        nonreq_columns = self.data.columns.drop(REQ_COLUMNS).tolist()
        output_columns = list(REQ_COLUMNS) + (nonreq_columns)
        self.data = self.data[output_columns]


    def process_file(self):
        """A set of main functions in Infiler to process input files.

        Following actions are performed:

        01. Check if file exists in the indicated download directory
        02. Load
        For each chunk,
          03. Check if empty.
          04. Standardize column names.
          05. Standardize chromosome names.
          06. Check columns for their expected variable data types.
          07. If file data_type == eqtl, check feature_id column.
             Standardize feature_id.
          08. Drop duplicate columns, if any.
          09. Reset index and reorder mandatory columns to the front.
          10. Write to file.
        """
        self.file_name = find_file(self.file_name)

        logging.info(f' loading {self.file_name}.')
        chunks = pandas.read_csv(
            self.file_name,
            # c engine does not support regex
            # sep=r'\s{,8}',
            # lots of files are whitespace delimited
            # default behavior will push all missing columns to last
            # delim_whitespace=True,
            # sep='\t| ',
            sep='\t',
            header=0,
            engine='python',
            iterator=True,
            index_col=None,
            chunksize=self.chunksize
        )

        chunkcount = 0

        for chunk in chunks:
            logging.debug(f' processing data.head(2): {chunk.head(2)}.')
            chunk.reset_index(drop=True, inplace=True)
            logging.info(f' processing input chunk {chunkcount}.')

            # check if empty and check header
            if not chunk.empty:
                if self.columnset:
                    self.rename_columns(chunk)
                    if '#CHROM' in chunk.columns:
                        logging.debug(f' renaming #CHROM to variant_chrom.')
                        chunk.rename(
                            columns={'#CHROM':'variant_chrom'},
                            inplace=True
                        )
                # check each column for variable types,
                # standardize chromosome and variant ids, etc.
                self.check_data(chunk)
                logging.debug(f' processing data type {self.data_type}.')

                if self.data_type == 'eqtl':
                    features = self.list_features()
                    # 413 error
                    # self.call_querier(features)
                    self.map_features(features)

                logging.info(f' dropping duplicate columns.')
                dropcols = self.data.columns.duplicated()
                self.data = self.data.loc[:, ~dropcols]

                # reorder columns so that mandatory fields are listed first.
                logging.info(f' reordering processed data.')
                self.order_columns()
                logging.debug(f' data.head(2): {self.data.head(2)}')

                logging.info(f' writing processed data.')

                # write if first chunk. append if not.
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


