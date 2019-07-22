"""Set of default values used in cimr"""


DATA_TYPES = ('gwas', 'eqtl', 'sqtl', 'pqtl', 'gene', 'twas', 'tad')

GENOME_BUILDS = ('b37', 'b38')

WORKING_HEADERS = ['rsnum', 'variant_id', 'pvalue', 
    'effect_size', 'odds_ratio', 'standard_error', 'zscore', 
    'tss_distance', 'effect_allele', 'non_effect_allele', 
    'frequency', 'imputation_status', 'sample_size', 'n_cases',
    'chrom', 'pos', 'ref', 'alt', 'build', 'chromosome', 
    'position', 'inc_allele', 'feature_id', 'feature_start',
    'feature_stop', 'feature_type', 'feature_name'
]

HEADERS = ['rsnum', 'variant_id', 'pvalue', 'effect_size', 
    'odds_ratio', 'standard_error', 'zscore', 'tss_distance', 
    'effect_allele', 'non_effect_allele', 'frequency', 
    'imputation_status', 'sample_size', 'n_cases'
]

