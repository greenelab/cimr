<<<<<<< HEAD
"""Set of default values used in cimr"""

=======
'''Set of default values used in cimr'''
>>>>>>> 33f450c038214b0e060f5f0899d82fb8f8107df4

DATA_TYPES = ('gwas', 'eqtl', 'sqtl', 'pqtl', 'gene', 'twas', 'tad')

GENOME_BUILDS = ('b37', 'b38')

<<<<<<< HEAD
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

=======
HEADERS = ['gene_id', 'rsnum', 'variant_id', 'pvalue', 
    'effect_size', 'standard_error', 'zscore', 'tss_distance', 
    'ma_samples', 'maf', 'chrom', 'pos', 'ref', 'alt', 
    'build', 'effect_allele', 'non_effect_allele', 'inc_allele', 
    'inc_afrq', 'imputation_status', 'sample_size', 'n_cases', 
    'frequency'
]
>>>>>>> 33f450c038214b0e060f5f0899d82fb8f8107df4
