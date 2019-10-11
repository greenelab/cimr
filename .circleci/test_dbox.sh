#!/bin/bash

set -ex

cd ~/cimr/examples/dbox/

# Run example yaml file for a dummy data
cimr processor -process -yaml-file dbox_gwas.yml

# Decompress for comparisons with arbitrary file change
# due to gzip defaults
gunzip -k processed_data/gwas/*.gz

# Verify md5 hash for processed file
md5sum --check md5_dbox_gwas.txt

# Check catalog.txt
diff catalog.txt expected_catalog.txt
rm -f catalog.txt

