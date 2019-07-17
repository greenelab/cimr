'''Set of default values used in cimr'''

DATA_TYPES = ('gwas', 'eqtl', 'sqtl', 'pqtl', 'gene', 'twas', 'tad')

GENOME_BUILDS = ('b37', 'b38')

HEADERS = ['gene_id', 'rsnum', 'variant_id', 'pvalue', 
    'effect_size', 'standard_error', 'zscore', 'tss_distance', 
    'ma_samples', 'maf', 'chrom', 'pos', 'ref', 'alt', 
    'build', 'effect_allele', 'non_effect_allele', 'inc_allele', 
    'inc_afrq', 'imputation_status', 'sample_size', 'n_cases', 
    'frequency'
]