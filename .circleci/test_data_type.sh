#!/bin/bash

set -e -x

cd ~/cimr/examples/

# Process yaml examples
for DATATYPE in "eqtl" "gwas"; do
    # Get in the specific dir
    cd ${DATATYPE}/

    # Process example yaml file
    cimr processor -process -yaml-file ${DATATYPE}.yml

    OUTPUT_PREFIX=${DATATYPE}/${DATATYPE}
    gunzip -k processed_data/${OUTPUT_PREFIX}.txt.gz

    # Verify output result
    md5sum --check md5_${DATATYPE}.txt

    # Delete "submitted_data" and "processed_data" directories
    rm -rf submitted_data processed_data

    # Verify catalog file
    diff catalog.txt expected_catalog.txt
    rm -f catalog.txt

    # Go back to examples dir
    cd ..
done


