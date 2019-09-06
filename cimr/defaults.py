"""Set of default values used in cimr"""

ANNOTURL = 'https://raw.githubusercontent.com/greenelab/cimr/master/cimr/data/annotation/'

MAXCHROM = 27
CHUNKSIZE = 7000000

VAR_COMPONENTS = [
    'variant_chrom',
    'variant_pos',
    'non_effect_allele',
    'effect_allele'
]

DATA_TYPES = (
    'gwas', 'eqtl', 'sqtl', 'pqtl', 'gene', 'twas',
    'tad', 'multiple', 'yaml'
)

GENOME_BUILDS = ('b37', 'b38')

WORKING_HEADER = [
    'rsnum', 'variant_id', 'pvalue',
    'effect_size', 'odds_ratio', 'standard_error', 'zscore',
    'tss_distance', 'effect_allele', 'non_effect_allele',
    'frequency', 'imputation_status', 'sample_size', 'n_cases',
    'chrom', 'pos', 'ref', 'alt', 'build', 'chromosome',
    'position', 'inc_allele', 'feature_id', 'feature_start',
    'feature_stop', 'feature_type', 'feature_name'
]

# required header
REQ_HEADER = [
    'variant_id', 'pvalue', 'effect_size', 'standard_error'
]

# union of expected columns in variant- or variant-feature
# association test results
HEADER = [
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
]

CONFIG_FILE_EXTENSION = ('yml', 'yaml')
# transparent compressions recognized by tarfile 'r:*' are:
# gzip, bz2, and lzma (xz)
COMPRESSION_EXTENSION = ('gz', 'bz2', 'xz')
BULK_EXTENSION = ('tgz', 'tar.gz', 'tar.bz2', 'tar.xz')
FILE_EXTENSION = ('txt.gz', 'tsv.gz', 'tsv', 'txt')

META_HEADER = [
    'data_type', 'context', 'context_id',
    'file_name', 'processed_data_url', 'submitted_data_url',
    'submitted_data_md5', 'citation', 'data_source',
    'build', 'sample_size', 'n_cases', 'method_name', 'method_tool',
    'description'
]
