#!/bin/bash

# Remove previous output files
rm -rf processed_data/ catalog.txt foo_*

# cimr processing
# default chunk size: 5 million
#cimr processor -process -yaml-file gwas.yml

# Set chunk size to 1 million
cimr processor -chunksize 1000000 -process -yaml-file gwas.yml

if [ -f processed_data/output.txt.gz ]; then        # serial processing
    gunzip processed_data/output.txt.gz
else                                            # parallel processing
    echo "`date`: combining output files ..."
    cat foo_* >> processed_data/output.txt
    echo "`date`: gzip output file ..."
    gzip --keep processed_data/output.txt
fi

# md5sum check
echo "`date`: md5sum check ..."
md5sum --check md5.txt

echo "`date`: DONE"
