"""Set of default values used in cimr"""

MAXCHROM = 23

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

HEADER = ['rsnum', 'variant_id', 'pvalue', 'effect_size', 
    'odds_ratio', 'standard_error', 'zscore', 'tss_distance', 
    'effect_allele', 'non_effect_allele', 'frequency', 
    'imputation_status', 'sample_size', 'n_cases', 'build'
]

CONFIG_FILE_EXTENSION = ('yml', 'yaml')
# transparent compressions recognized by tarfile 'r:*' are:
# gzip, bz2, and lzma (xz)
COMPRESSION_EXTENSION = ('gz', 'bz2', 'xz')
BULK_EXTENSION = ('tgz', 'tar.gz', 'tar.bz2', 'tar.xz')
FILE_EXTENSION = ('txt.gz', 'tsv.gz')

META_HEADER = [
    'file_name', 'processed_data_url', 'submitted_data_url', 
    'submitted_data_md5', 'citation', 'data_source', 'data_type', 'build',
    'sample_size', 'n_cases', 'method_name', 'method_tool', 'description'
]
