#!/bin/bash
#
# This script tests gwas data processing with different "chunksize" and
# "parallel" parameters.
#
# Note that in this test, the decompressed input GWAS file has 2,201 rows
# (header line included).

set -ex

DIRNAME=`dirname $0`
cd $DIRNAME
rm -rf submitted_data/ processed_data/ catalog.txt

echo "Test with default chunksize and parallel parameters"
cimr processor -process -yaml-file gwas.yml
gunzip -k processed_data/gwas/gwas.txt.gz
diff catalog.txt expected_output/catalog.txt
diff processed_data/gwas/gwas.txt expected_output/gwas.txt
rm -rf processed_data/ catalog.txt
echo

echo "Test with chunksize=500 and parallel=2"
cimr processor -process -yaml-file gwas.yml -chunksize 500 -parallel 2
gunzip -k processed_data/gwas/gwas.txt.gz
diff catalog.txt expected_output/catalog.txt
diff processed_data/gwas/gwas.txt expected_output/gwas.txt
rm -rf processed_data/ catalog.txt
echo

echo "Test with default chunksize and non-parallel (parallel=0)"
cimr processor -process -yaml-file gwas.yml -parallel 0
gunzip -k processed_data/gwas/gwas.txt.gz
diff catalog.txt expected_output/catalog.txt
diff processed_data/gwas/gwas.txt expected_output/gwas.txt
rm -rf processed_data/ catalog.txt
echo

echo "Test with chunksize=500 and non-parallel (parallel=0)"
cimr processor -process -yaml-file gwas.yml -chunksize 500 -parallel 0
gunzip -k processed_data/gwas/gwas.txt.gz
diff catalog.txt expected_output/catalog.txt
diff processed_data/gwas/gwas.txt expected_output/gwas.txt
rm -rf processed_data/ catalog.txt
echo

# Remove downloaded file
rm -rf submitted_data
echo "Tests successful"
