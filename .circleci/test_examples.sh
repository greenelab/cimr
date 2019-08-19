#!/bin/bash

set -e -x

cd examples/

# Process yaml examples
for DATATYPE in "eqtl" "gwas"; do
    cimr processor -process -yaml-file submitted/${DATATYPE}.yml

    OUTPUT_PREFIX=${DATATYPE}/${DATATYPE}
    gunzip -k processed_data/${OUTPUT_PREFIX}.txt.gz

    # Verify output result
    diff processed_data/${OUTPUT_PREFIX}.txt expected_processed_data/${OUTPUT_PREFIX}.txt

    # Clear up "submitted_data" and "processed_data"
    rm -rf submitted_data processed_data
done

# Verify catalog file
diff cimr-d_catalog.txt expected_cimr-d_catalog.txt
