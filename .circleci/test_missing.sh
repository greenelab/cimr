#!/bin/bash

set -e -x

cd ~/cimr/examples/

# Tested data type
DATA_TYPE='gwas'

# Process yaml examples
for MISSING in "standard_error" "pvalue" "zscore"; do
    # Test data with missing data
    cd ${MISSING}/

    cimr processor -process -yaml-file ${DATA_TYPE}_no_${MISSING}.yml
    OUTPUT_PREFIX=${DATA_TYPE}/${DATA_TYPE}_no_${MISSING}
    gunzip -k processed_data/${OUTPUT_PREFIX}.txt.gz

    # Verify output result
    md5sum --check md5_${DATATYPE}_no_${MISSING}.txt

    # Delete "submitted_data" and "processed_data" directories
    rm -rf submitted_data processed_data

    # Verify catalog file
    diff catalog.txt expected_catalog.txt
    rm -f catalog.txt

    # Go back to examples dir
    cd ..
done


