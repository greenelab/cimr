


cimr processor -process \
    -file-name example_submitted/gwas.txt.gz \
    -data-type gwas \
    -column-set '{"variant_id":"panel_variant_id", "variant_chrom":"chromosome", "variant_pos":"position", "rsnum":"variant_id"}' \
    -outdir example_processed \
    -out gwas


cimr processor -process \
    -file-name example_submitted/eqtl.txt.gz \
    -data-type eqtl \
    -column-set '{"effect_size":"slope", "standard_error":"slope_se", "pvalue":"pval_nominal", "feature_id":"gene_id"}' \
    -outdir example_processed \
    -out eqtl


cimr processor -process \
    -file-name example_submitted/sqtl.txt.gz \
    -data-type sqtl \
    -column-set '{"rsnum":"rs_id_dbSNP151_GRCh38p7", "effect_size":"slope", "standard_error":"slope_se", "pvalue":"pval_nominal", "feature_id":"gene_id", "feature_start":"gene_start", "feature_stop":"gene_end", "feature_name":"gene_name", "feature_chrom":"gene_chr", "chrom":"chr"}' \
    -outdir example_processed \
    -out sqtl



