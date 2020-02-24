#!/usr/bin/env python3

"""cimr processor submodule main classes

(c) YoSon Park
"""

import copy
import multiprocessing as mp
import os
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
    check_numeric, check_probability, remove_palindromic)

# default values
from ..defaults import (COMPRESSION_EXTENSION, ANNOT_DIR,
    DATA_TYPES, GENOME_BUILDS, VAR_COMPONENTS, SEPARATORS,
    MAXCHROM, HEADER, REQ_COLUMNS, NUMERIC_COLUMNS,
    PROB_COLUMNS, INT_COLUMNS,
    SNP125HG17, SNP130HG18, SNP150HG19, SNP150HG38, HG19TO38)


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
                 parallel,
                 columnset):

        if data_type not in DATA_TYPES:
            raise ValueError(
                ' %s is not a valid data_type supported' % data_type
            )
        if genome_build not in GENOME_BUILDS:
            raise ValueError(
                ' %s is not a valid genome_build supported' % genome_build
            )
        self.data = pandas.DataFrame()
        self.data_type = data_type
        self.file_name = file_name
        self.genome_build = genome_build
        self.update_rsid = update_rsid
        self.outfile = outfile
        self.chunksize = chunksize
        self.parallel = parallel
        self.columnset = columnset


    def get_sep(self):
        """Get column separator, if provided."""
        if 'column_separator' in self.columnset.values():
            logging.debug(f' column_separator is provided.')
            _keys = list(self.columnset.keys())
            _values = list(self.columnset.values())
            sep = _keys[_values.index('column_separator')]
            sep = sep.upper()

            if sep in SEPARATORS:
                self.sep = SEPARATORS[sep]
            else:
                logging.warning(f' using tab as the default.')
                self.sep = '\t'

        else:
            self.sep = '\t'


    def get_pos(self):
        """Check variant_id column and make
        chrom, pos, ref, alt, build columns"""
        # check first element to find delimiter in variant_id
        variant_ids = self.data['variant_id']
        if '_' in variant_ids[0]:
            pass
        if '-' in variant_ids[0]:
            variant_ids = variant_ids.str.replace('-', '_')
            logging.info(
                f' chunk #{self.chunk_id}:' +
                f' variant_id delimiter changed from \'-\' to \'_\''
            )
        if ':' in variant_ids[0]:
            variant_ids = variant_ids.str.replace(':', '_')
            logging.info(
                f' chunk #{self.chunk_id}:' +
                f' variant_id delimiter changed from \':\' to \'_\''
            )
        # else:
        #     logging.error(
        #        f' unknown delimiter used in variant_id'
        #     )
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
                logging.info(
                    f' chunk #{self.chunk_id}:' +
                    f' updating variant_id to include build'
                )
                self.data['build'] = self.genome_build
                self.data['variant_id'] = variant_ids + '_' + self.genome_build


    def make_variant_id(self):
        """(Re)make variant_id with updated chrom ID, build #, etc."""
        logging.debug(
            f' chunk #{self.chunk_id}: standardizing allele codes...'
        )
        self.data['ref'] = self.data['ref'].str.upper()
        self.data['alt'] = self.data['alt'].str.upper()

        logging.debug(
            f' chunk #{self.chunk_id}: making a new variant_id column.'
        )

        self.data['variant_id'] = self.data['chrom'].astype(str) \
            + '_' \
            + self.data['pos'].astype(str) \
            + '_' \
            + self.data['ref'].astype(str) \
            + '_' \
            + self.data['alt'].astype(str) \
            + '_' \
            + self.data['build'].astype(str)
        logging.info(
            f' chunk #{self.chunk_id}: variant_id has been standardized.'
        )


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
            logging.debug(f' chunk #{self.chunk_id}: {chroms}')

        if len(chroms) > (maxchrom - 2) and len(chroms) < (maxchrom + 2):
            logging.info(
                f' chunk #{self.chunk_id}: there are {len(chroms)} chromosomes.'
            )
        elif len(chroms) <= (maxchrom - 2):
            logging.info(
                f' chunk #{self.chunk_id}:' +
                f' chromosome(s) included: %s' % (chroms,)
            )
        else:
            logging.warning(
                f' chunk #{self.chunk_id}:' +
                f' input file more than {maxchrom - 1} chromosomes.'
            )
            logging.warning(
                f' chunk #{self.chunk_id}: chromosome(s) included: %s'%(chroms,)
            )

        if len(set(chroms) & set(chrom_str)) > (maxchrom - 2):
            pass

        if len(set(chroms) & set(chrom_int)) > 0:
            self.data['chrom'] = self.data['chrom'].map(
                chrom_dict, na_action='ignore'
            ).fillna(self.data['chrom'])
            logging.info(
                f' chunk #{self.chunk_id}: chromosome ids have been updated.'
            )
        else:
            logging.warning(
                f' chunk #{self.chunk_id}: chromosome id needs to be checked.'
            )

        remainder = list(set(chroms) - set(chrom_str) - set(chrom_int))
        if len(remainder) > 0:
            logging.warning(
                f' chunk #{self.chunk_id}:' +
                f' chromosome(s) not used: %s' % (remainder,)
            )

        logging.info(
            f' chunk #{self.chunk_id}: finished checking chromosome ids.'
        )


    def check_ref(self):
        """Compare a random subset of variants to GRCh37 or GRCh38 references."""
        from pkg_resources import resource_filename

        var_ref = self.var_ref
        ref_id = self.var_ref_id
        ref_file = resource_filename(
            'cimr', var_ref
            )
        logging.info(
            f' chunk #{self.chunk_id}: using {ref_file} to check variants.'
        )
        logging.info(f' chunk #{self.chunk_id}: rs id reference is {ref_id}.')
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
            variant_nomatch = merged.loc[
                ~(merged['variant_id']==merged['variant_id_reference'])
            ]
            samplecount = len(samples.index)
            rsrefcount = len(rsnum_ref.index)
            nomatchcount = len(variant_nomatch.index)
            logging.info(
                f' chunk #{self.chunk_id}:' +
                f' {samplecount} sampled from {rsrefcount} total' +
                f' variants with rs ids,'
            )
            logging.info(
                f' chunk #{self.chunk_id}:' +
                f' {nomatchcount} variants do not match the reference.'
            )
        else:
            logging.error(
                f' chunk #{self.chunk_id}: there are no matching rs ids.'
            )
            pass


    def find_reference(self):
        """Find variant and gene references for map checking"""
        if self.genome_build == 'b37':
            self.gene_ref = os.path.join(
                ANNOT_DIR, 'gene_grch37_gencode_v26.tsv.gz'
            )
            self.var_ref = os.path.join(
                ANNOT_DIR, 'variant_grch37_subset.txt.gz'
            )
            self.var_ref_id = 'rs_id_dbSNP147_GRCh37p13'
            self.dbsnp = SNP150HG19
            self.chain = HG19TO38
        elif self.genome_build == 'b38':
            self.gene_ref = os.path.join(
                ANNOT_DIR, 'gene_grch38_gencode_v29.tsv.gz'
            )
            self.var_ref = os.path.join(
                ANNOT_DIR, 'variant_grch38_subset.txt.gz'
            )
            self.var_ref_id = 'rs_id_dbSNP150_GRCh38p7'
            self.dbsnp = SNP150HG38
        else:
            logging.error(f' accepted genome_build values: b37, b38.')
            sys.exit(1)


    def trim_ensembl(self):
        """Make the ensembl gene id list to be agnostic to
        subversions for queries

        e.g. ENSG00000143614.7 -> ENSG00000143614
        """
        self.data['original_feature_id'] = self.data['feature_id']
        ensemblid = self.data['feature_id'].str.split(r'\.').str[0]
        self.data['feature_id'] = ensemblid
        logging.info(
            f' chunk #{self.chunk_id}:' +
            f' ensembl id has been truncated for database queries.'
        )


    def list_features(self):
        """Find the list of features (e.g. genes)."""
        logging.debug(f' chunk #{self.chunk_id}: {self.included_header}')
        logging.debug(f' chunk #{self.chunk_id}: {self.data.head(2)}')
        if 'feature_id' in self.included_header:
            if self.data['feature_id'].iloc[0].startswith('ENS'):
                self.trim_ensembl()

        if 'ensemblgene' in self.data.columns:
            return self.data.ensemblgene
        elif 'feature_id' in self.data.columns:
            return self.data.feature_id
        else:
            logging.error(
                f' chunk #{self.chunk_id}: feature_id column is not provided.'
            )
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
        logging.debug(
            f' chunk #{self.chunk_id}:' +
            f' processing information needed to make variant_id.'
        )

        if set(VAR_COMPONENTS).issubset(set(self.included_header)):
            logging.debug(
                f' chunk #{self.chunk_id}:' +
                f' checking {VAR_COMPONENTS} for missing values.'
            )
            self.data.dropna(subset=VAR_COMPONENTS, inplace=True)
            logging.debug(
                f' chunk #{self.chunk_id}:' +
                f' rows with missing {VAR_COMPONENTS} are dropped.'
            )

            logging.debug(
                f' chunk #{self.chunk_id}:' +
                f' making variant_id using {VAR_COMPONENTS}.'
            )
            self.data['chrom'] = self.data['variant_chrom']
            self.check_chrom()
            logging.info(
                f' chunk #{self.chunk_id}:' +
                f' verifying variant positions are int values.'
            )
            self.data['variant_pos'] = self.data['variant_pos'].astype(int)
            logging.debug(
                f' chunk #{self.chunk_id}:' +
                f' changing verified values to str.'
            )
            self.data['variant_pos'] = self.data['variant_pos'].astype(str)

            self.data['pos'] = self.data['variant_pos']
            self.data['ref'] = self.data['non_effect_allele']
            self.data['alt'] = self.data['effect_allele']
            self.data['build'] = self.genome_build

            # make variant_id from components
            self.make_variant_id()

            # redefine included_header from columns including variant_id
            self.included_header = intersect_set(
                HEADER, self.data.columns
            )

        else:
            logging.error(
                f' chunk #{self.chunk_id}:' +
                f' missing columns necessary to make variant_id.'
            )
            sys.exit(1)


    def check_data(self, chunk):
        """Check different columns for dtype, remove missing rows and
        standardize format to be used for analyses
        """
        logging.debug(f' chunk #{self.chunk_id}: HEADER: {HEADER}.')
        logging.debug(f' chunk #{self.chunk_id}: columns: {chunk.columns}.')
        self.included_header = intersect_set(
            HEADER, chunk.columns
        )
        logging.debug(
            f' chunk #{self.chunk_id}:' +
            f' included header overlapping cimr default set:' +
            f' {self.included_header}.'
        )

        self.data = chunk.copy()  # self.data is a deep copy of chunk
        self.data.reset_index(inplace=True, drop=True)

        # standardize variant_id based on info provided in the input file
        if 'variant_id' not in self.data.columns:
            logging.info(
                f' chunk #{self.chunk_id}: variant_id column is not provided.'
            )
            logging.info(
                f' chunk #{self.chunk_id}: ' +
                'checking columns necessary to make variant_id...'
            )
            self.make_composite_id()
            logging.debug(
                f' chunk #{self.chunk_id}: chunk.head(2): {self.data.head(2)}'
            )
        elif 'variant_id' in self.included_header:
            self.get_pos()
            self.check_chrom()
            self.make_variant_id()
        else:
            logging.error(
                f' chunk #{self.chunk_id}: variant_id column is not provided.'
            )
            sys.exit(1)

        # update map + variant_id if genome_build is not hg38
        if (self.genome_build != 'hg38') and (self.genome_build != 'b38'):
            from .lift import call_liftover
            self.data = call_liftover(self.data)
            self.make_variant_id()

        if 'rsnum' in self.data.columns:
            if self.update_rsid:
                self.check_ref()
        else:
            logging.warning(
                f' chunk #{self.chunk_id}: rsnum column is not provided.'
            )

        for col in INT_COLUMNS:
            if col in self.data.columns:
                self.data[col] = self.data[col].astype(int)

        if 'inc_allele' in self.data.columns:
            self.fill_effect_allele()
        else:
            logging.debug(
                f' chunk #{self.chunk_id}: inc_allele column is not available.'
            )


        self.data = remove_palindromic(self.data)
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
                logging.error(f' chunk #{self.chunk_id}: {col} is required.')
                sys.exit(1)

        for col in NUMERIC_COLUMNS:
            if col in self.data.columns:
                self.data = check_numeric(self.data, col)

        for col in PROB_COLUMNS:
            if col in self.data.columns:
                check_probability(self.data, col)


    def process_chunk(self, chunk, err_queue=None):
        try:
            logging.info(f' chunk #{self.chunk_id}: start processing ...')
            logging.debug(
                f' chunk #{self.chunk_id}:' +
                f' processing data.head(2): {chunk.head(2)}.'
            )
            chunk.reset_index(drop=True, inplace=True)

            # check if empty and check header
            if not chunk.empty:
                if self.columnset:
                    self.rename_columns(chunk)
                    if '#CHROM' in chunk.columns:
                        logging.debug(
                            f' chunk #{self.chunk_id}:' +
                            f' renaming #CHROM to variant_chrom.'
                        )
                        chunk.rename(
                            columns={'#CHROM':'variant_chrom'},
                            inplace=True
                        )
                # check each column for variable types,
                # standardize chromosome and variant ids, etc.
                self.check_data(chunk)
                logging.debug(
                    f' chunk #{self.chunk_id}:' +
                    f' processing data type {self.data_type}.'
                )

                if self.data_type == 'eqtl':
                    features = self.list_features()
                    # 413 error
                    # self.call_querier(features)
                    self.map_features(features)

                logging.info(
                    f' chunk #{self.chunk_id}: dropping duplicate columns.'
                )
                dropcols = self.data.columns.duplicated()
                self.data = self.data.loc[:, ~dropcols]

                # reorder columns so that mandatory fields are listed first.
                logging.info(
                    f' chunk #{self.chunk_id}: reordering processed data.'
                )
                self.order_columns()
                logging.debug(
                    f' chunk #{self.chunk_id}:'
                    f' data.head(2): {self.data.head(2)}'
                )

                logging.info(
                    f' chunk #{self.chunk_id}: writing processed data.'
                )
                self.write_chunk()
                logging.info(
                    f' chunk #{self.chunk_id}: ' +
                    '-' * 22 + ' DONE ' + '-' * 22
                )
            else:
                logging.error(
                    f' chunk #{self.chunk_id}: no content in {self.file_name}.'
                )
                sys.exit(1)
        except Exception as e:
            err_message = f' chunk #{self.chunk_id}: error'
            if self.parallel > 0:
                import traceback
                err_queue.put(err_message + '\n' + traceback.format_exc())
            else:
                # In non-parallel mode, log the error message, and propagate
                # the exception to upper level of the program.
                logging.error(err_message)
                raise


    def write_chunk(self):
        """Write chunk data to a chunk file without header"""

        is_chunk1 = (self.chunk_id == 1)
        if is_chunk1 or self.parallel > 0:
            mode = 'w'
        else:
            mode = 'a'

        if self.parallel > 0:
            out_filename = self.chunk_file_prefix + str(self.chunk_id)
        else:
            out_filename = str(self.outfile)

        if isinstance(self.data, pandas.DataFrame):
            self.data.to_csv(
                out_filename,
                header=is_chunk1,
                index=False,
                sep='\t',
                na_rep='NA',
                compression='gzip',
                float_format='%.6f',
                mode=mode,
                encoding='utf-8'
            )


    def combine_chunks(self):
        with open(str(self.outfile), 'wb') as out_file:
            for i in range(1, self.num_chunks + 1):
                chunk_filename = self.chunk_file_prefix + str(i)
                with open(chunk_filename, 'rb') as in_file:
                    out_file.write(in_file.read())
                os.remove(chunk_filename)  # remove chunk output file

    def check_h5_outfile(self):
        """Check file extension for h5 output"""
        if not str(self.outfile).endswith('.h5.bz2'):
            self.outfile = str(self.outfile) + '.h5.bz2'

    def rename_columns(self, chunk):
        """Given a pandas dataframe and a dictionary, rename columns
        Primarily used for yaml column-set-based renaming in cimr.
        """
        logging.info(
            f' chunk #{self.chunk_id}:' +
            f' renaming columns based on column dictionary.'
        )
        chunk.rename(self.columnset, axis=1, inplace=True)
        logging.debug(
            f' chunk #{self.chunk_id}: renamed data.head(2): {chunk.head(2)}'
        )


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
            logging.error(
                f' chunk #{self.chunk_id}:' +
                f' feature reference cannot be determined.'
            )
            sys.exit(1)

        logging.debug(f' chunk #{self.chunk_id}: reading {self.gene_ref}.')

        if not 'feature_name' in self.included_header:
            gene_annot = pandas.read_csv(
                self.gene_ref,
                sep='\t',
                header=0
            )
            cols = ['feature_name', self.feature_type, 'feature_type']
            gene_annot = gene_annot[cols]
            gene_annot.rename(
                columns={self.feature_type:'feature_id'}, inplace=True
            )

            logging.debug(
                f' chunk #{self.chunk_id}: ' +
                f'current dataframe: {self.data.head(2)}'
            )
            logging.debug(
                f' chunk #{self.chunk_id}:' +
                f' selected annotations: {gene_annot.head(2)}'
            )

            self.data = self.data.merge(
                gene_annot,
                on='feature_id',
                how='left',
                left_index=False,
                right_index=False
            )
            logging.debug(
                f' chunk #{self.chunk_id}: dropping duplicates, if any...'
            )
            self.data.drop_duplicates(inplace=True)
            self.data.reset_index(inplace=True, drop=True)
            logging.debug(
                f' chunk #{self.chunk_id}: dataframe has been annotated.'
            )
            logging.debug(f' chunk #{self.chunk_id}: {self.data.head(2)}')


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
          06. If variant_id is not provided or non-standard format,
              search for chromosome, position, ref, alt and genome_build
              information to create a standardized variant_id column
          07. If effect_size, pvalue or zscore columns are not provided,
              search for necessary information to estimate using
              convertibles
          08. Check columns for their expected variable data types.
          09. If file data_type == eqtl, check feature_id column.
              Standardize feature_id.
          10. Drop duplicate columns, if any.
          11. Reset index and reorder mandatory columns to the front.
          12. Write to file.
        """

        self.file_name = find_file(self.file_name)
        self.find_reference()

        logging.debug(f' looking for column separators.')
        self.get_sep()

        logging.info(f' loading {self.file_name}.')
        chunks = pandas.read_csv(
            self.file_name,
            # c engine does not support regex
            # sep=r'\s{,8}',
            # lots of files are whitespace delimited
            # default behavior will push all missing columns to last
            # delim_whitespace=True,
            # sep='\t| ',
            sep=self.sep,
            header=0,
            iterator=True,
            index_col=None,
            chunksize=self.chunksize
        )

        logging.info(f' chunksize: {self.chunksize / 1000000} million')
        logging.info(f' number of parallel processes: {self.parallel}')

        chunkcount = 0
        if self.parallel > 0:  # parallel data processing
            self.chunk_file_prefix = str(self.outfile) + ".chunk"
            pool = mp.Pool(self.parallel)
            err_queue = mp.Manager().Queue()  # a queue shared by all processes
            for chunk in chunks:
                if not err_queue.empty():
                    break    # break immediately if err_queue is not empty
                chunkcount += 1
                self.chunk_id = chunkcount
                cloned_instance = copy.deepcopy(self)
                pool.apply_async(
                    cloned_instance.process_chunk,
                    [chunk, err_queue]
                )
            pool.close()
            pool.join()

            # Terminate the whole program if error is found in child processes.
            # Because child processes are launched in parallel, err_queue may
            # have multiple entries, but we only care the first error.
            if not err_queue.empty():
                logging.error(err_queue.get())
                sys.exit(1)

            # Combine all chunk output files into a single output file
            logging.info(f" combine chunk output files ...")
            self.num_chunks = chunkcount
            self.combine_chunks()

        else:  # serial (non-parallel) data processing
            chunkcount = 0
            for chunk in chunks:
                chunkcount += 1
                self.chunk_id = chunkcount
                self.process_chunk(chunk)


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
            clonercmd = (
                'git clone git@github.com:greenelab/cimr-d.git ' + self.tempdir
            )
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
