"""Set of default values used in cimr"""


import os
from collections import namedtuple

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ANNOT_DIR = os.path.join(BASE_DIR, 'data', 'annotation')
CHAIN_DIR = os.path.join(BASE_DIR, 'data', 'chain')

# Current maximum is set with human chromosomes:
# 1 - 22 autosomal chromosomes

# plink and related software treat other chromosomes as:
# X chromosome (n+1; 23)
# Y chromosome (n+2; 24)
# XY chromosome (n+3; 25)
# MT chromosome (n+4; 26)
MAXCHROM = 27

# Chunk size to read in submitted data
# Current default is set to accomodate complex format data
# (with dtype=object columns in pandas)
# in CI environment with limited memory
# Can be overwritten with commandline argument '-chunksize'
CHUNKSIZE = 2000000

# Default number of parallel processes,
# can be overwritten with commandline argument '-parallel'
PARALLEL = 4

VERY_SMALL_P = 1e-70

# Minimally sufficient information to make unique variant IDs
VAR_COMPONENTS = [
    'variant_chrom',
    'variant_pos',
    'non_effect_allele',
    'effect_allele'
]

# Column separators commonly found in input files
SEPARATORS = {
    'TAB': '\t',
    'SPACE': ' ',
    'TABSPACE': '\t| ',
    'MULTISPACE': r'\s+',
    '8SPACE': r'\s{,8}',
    'COMMA': ','
}

DATA_TYPES = (
    'gwas', 'eqtl', 'sqtl', 'pqtl', 'gene', 'twas',
    'tad', 'multiple', 'yaml'
)

GENOME_BUILDS = ('b37', 'b38')

# temporary local file locations
SNP125HG17 = '/data/refs/dbsnp/snp125_hg17.txt.gz'
SNP130HG18 = '/data/refs/dbsnp/snp130_hg18.txt.gz'
SNP150HG19 = '/data/refs/dbsnp/snp150_hg19.txt.gz'
SNP150HG38 = '/data/refs/dbsnp/snp150_hg38.txt.gz'

HG19TO38 = os.path.join(CHAIN_DIR, 'hg19ToHg38.over.chain.gz')
HG18TO38 = os.path.join(CHAIN_DIR, 'hg18ToHg38.over.chain.gz')

DBSNP_HEADER = namedtuple('DBSNP',
  ['binary_index', 'chromosome', 'start', 'end', 'name', 'score',
  'strand', 'ref_ncbi', 'ref_ucsc', 'observed', 'mol_type',
  'variant_class', 'valid', 'average_heterozigosity', 'average_heterozigosity_se',
  'functional_category', 'location_type', 'quality_weight', 'exceptions',
  'submitter_count', 'submitters', 'alleles_with_freq', 'alleles',
  'allele_chromosome_count', 'allele_frequencies', 'bitfields']
)

WORKING_HEADER = {
    'rsnum', 'variant_id', 'pvalue',
    'effect_size', 'odds_ratio', 'standard_error', 'zscore',
    'tss_distance', 'effect_allele', 'non_effect_allele',
    'frequency', 'imputation_status', 'sample_size', 'n_cases',
    'chrom', 'pos', 'ref', 'alt', 'build', 'chromosome',
    'position', 'inc_allele', 'feature_id', 'feature_start',
    'feature_stop', 'feature_type', 'feature_name'
}

# Required header
REQ_COLUMNS = [
    'variant_id', 'pvalue', 'effect_size', 'standard_error'
]

# Union of expected columns in variant- or variant-feature
# association test results
HEADER = {
    'variant_id', 'pvalue', 'effect_size', 'standard_error',
    'feature_id',
    'zscore', 'odds_ratio', 'rsnum', 'tss_distance',
    'effect_allele', 'non_effect_allele', 'frequency',
    'imputation_status', 'sample_size', 'n_cases', 'build',
    'feature_variant_id', 'feature_name',
    'feature_start', 'feature_stop', 'feature_type', 'strand',
    'ma_samples', 'ma_count', 'maf',
    'pvalue_perm', 'pvalue_beta', 'fdr', 'qvalue',
    'variant_chrom', 'variant_pos'
}

NUMERIC_COLUMNS = {
    'effect_size', 'standard_error', 'pvalue', 'pvalue_perm',
    'fdr', 'qvalue', 'odds_ratio', 'zscore'
}

INT_COLUMNS = {
    'ma_samples', 'sample_size', 'ma_count'
}

PROB_COLUMNS = {'pvalue', 'pvalue_perm', 'fdr'}

CONFIG_FILE_EXTENSION = ('yml', 'yaml')

# Transparent compressions recognized by tarfile 'r:*' are:
# gzip, bz2, and lzma (xz)
COMPRESSION_EXTENSION = ('gz', 'bz2', 'xz')
BULK_EXTENSION = ('tgz', 'tar.gz', 'tar.bz2', 'tar.xz')
FILE_EXTENSION = ('txt.gz', 'tsv.gz', 'tsv', 'txt')

META_HEADER = [
    'data_type', 'context', 'context_id',
    'file_name', 'processed_data_url', 'submitted_data_url',
    'submitted_data_md5', 'citation', 'data_source',
    'build', 'context_variable_type', 'sample_size', 'n_cases',
    'method_name', 'method_tool',
    'description', 'cimr_version'
]
