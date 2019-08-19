#!/bin/bash

set -e -x

cd examples/

# Process yaml examples
for DATATYPE in "eqtl" "gwas"; do
    cimr processor -process -yaml-file submitted/${DATATYPE}.yml

    OUTPUT_PREFIX=${DATATYPE}/${DATATYPE}
    gunzip -k processed_data/${OUTPUT_PREFIX}.txt.gz

    # Verify output result
    md5sum --check md5_${DATATYPE}.txt

    # Delete "submitted_data" and "processed_data" directories
    rm -rf submitted_data processed_data
done

# Verify catalog file
diff cimr-d_catalog.txt expected_cimr-d_catalog.txt
rm -f cimr-d_catalog.txt
