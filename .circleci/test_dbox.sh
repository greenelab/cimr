#!/bin/bash

set -ex

cd ~/cimr/examples/dbox/

# run example yaml file for a dummy data
cimr processor -process -yaml-file dbox_gwas.yml

# decompress for comparisons with arbitrary file change
# due to gzip defaults
gunzip -k processed_data/gwas/*.gz

# verify md5 hash for processed file
md5sum --check ../md5_dbox_gwas.txt
