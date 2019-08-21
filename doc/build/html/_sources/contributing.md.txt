
# How to contribute to cimr

cimr is an open source project that welcomes contributions from the
community, including:

* Data contributions

  * Variant association summary statistics
  * Annotations and references

* Method contributions

* Code contributions

  * Code patches,
  * Documentation improvements,
  * [Bug reports](https://github.com/greenelab/cimr/issues)
  * Reviews for contributed patches


## Did you find a bug?

Ensure the bug was not already reported by searching on
[github issues](https://github.com/greenelab/cimr/issues).
If there are no open issues addressing the problem(s), open a
new one. Be sure to clearly state a title, description and as much
relevant information as possible with example codes and logs.


## Did you fix a bug?

Open a new github pull request (PR) with the patch. Ensure that the
PR describes the problem and proposed solution. Refer to relevant
issues with a link. PR regarding cosmetic or stylistic edits are not
accepted.


## Do you have data?

Ensure that the data is publicly available or contains a license to
share with attribution (e.g. CC BY 4.0). Briefly, any variant- or
gene-based association summary statistics passing the automated checks
via cimr may be elligible to be integrated into cimr-d, which would
allow faster, more convenient use by others in the community.

All due credits, acknowledgements or citations are recommended for
users of this resource. cimr-d is provided as a public repository via
git. Information regarding versioned archives for data freezes will be
provided as they become available. License or other permissions may be
necessary to integrate into cimr-d.


### Variant-based association study results

Genome-wide association study (GWAS) and expression- or other
quantitative trait loci (eQTL, pQTL, sQTL, TWAS, etc.) have similar
input formatting requirements.


Following are required for all association summary statistics:

```
rsnum : rs id of the variant
variant_id : variant id in the following format
    chromosome_position_referenceallele_alternateallele_genomebuild
pval / pvalue: p-value of the beta estimate
beta / effect_size : beta coefficient estimate for the effect
    of the variant on the gene
se / standard_error : standard error of the beta
zscore: standardized or imputed statistic
sample_size : sample size for the variant underlying the
    effect_size and pvalue estimates. If not available per-variant,
    the total sample size used for the study is accepted.
n_cases : **required for binary phenotypic traits only.**
    n_cases are necessary
    for downstream analyses such as coloc colocalization.
build: if not included in variant_id, genome build can be provided
    as a yaml parameter
```

Following are required for quantitative trait loci:

```
feature_id : feature may be gene (eQTL), transcript (sQTL)
```

Following are optional:

```
tss_distance : distance to the transcription start site of the
    gene_id, if gene is used as a feature (eQTL)
ma_samples : samples with minor alleles
maf : minor allele frequency
effect_allele : allele with respect to which variant's effect
    sizes are estimated
non_effect_allele : non-effect-allele
imputation_status : imputation status of the variant. 'original'
    to indicate that the original GWAS contained the variant.
    'imputed' to indicate that the data submitter has done
    (additional) imputation to the initial release of the data.
frequency : allele frequency of the effect_allele
```

Following are assumed by parsing variant_id column:

```
chrom : chromosome id
pos: genomic position
ref : reference allele
alt : alternate allele
build : genomic build version number for the variant_id
```


