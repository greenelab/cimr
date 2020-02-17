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
    gunzip -dc processed_data/${OUTPUT_PREFIX}.txt.gz > processed_data/${OUTPUT_PREFIX}_default.txt

    cimr processor -process -yaml-file ${DATATYPE}.yml -chunksize 1000000
    gunzip -dc processed_data/${OUTPUT_PREFIX}.txt.gz > processed_data/${OUTPUT_PREFIX}_1m_chunksize.txt

    diff processed_data/${OUTPUT_PREFIX}_default.txt processed_data/${OUTPUT_PREFIX}_1m_chunksize.txt

    # Delete "submitted_data" and "processed_data" directories
    rm -rf submitted_data processed_data
    rm -f catalog.txt

    # Go back to examples dir
    cd ..
done

