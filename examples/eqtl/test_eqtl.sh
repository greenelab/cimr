#!/bin/bash
#
# This script tests gwas data processing with different "chunksize" and
# "parallel" parameters.
#
# Note that in this test, the decompressed input EQTL file has 40,119 rows
# (header line included).

set -ex

DIRNAME=`dirname $0`
cd $DIRNAME
rm -rf submitted_data/ processed_data/ catalog.txt

echo "Test with default chunksize and parallel parameters"
cimr processor -process -yaml-file eqtl.yml
gunzip -k processed_data/eqtl/eqtl.txt.gz
diff catalog.txt expected_output/catalog.txt
diff processed_data/eqtl/eqtl.txt expected_output/eqtl.txt
rm -rf processed_data/ catalog.txt
echo

echo "Test with chunksize=5000 and parallel=2"
cimr processor -process -yaml-file eqtl.yml -chunksize 5000 -parallel 2
gunzip -k processed_data/eqtl/eqtl.txt.gz
diff catalog.txt expected_output/catalog.txt
diff processed_data/eqtl/eqtl.txt expected_output/eqtl.txt
rm -rf processed_data/ catalog.txt
echo

echo "Test with default chunksize and non-parallel (parallel=0)"
cimr processor -process -yaml-file eqtl.yml -parallel 0
gunzip -k processed_data/eqtl/eqtl.txt.gz
diff catalog.txt expected_output/catalog.txt
diff processed_data/eqtl/eqtl.txt expected_output/eqtl.txt
rm -rf processed_data/ catalog.txt
echo

echo "Test with chunksize=5000 and non-parallel (parallel=0)"
cimr processor -process -yaml-file eqtl.yml -chunksize 5000 -parallel 0
gunzip -k processed_data/eqtl/eqtl.txt.gz
diff catalog.txt expected_output/catalog.txt
diff processed_data/eqtl/eqtl.txt expected_output/eqtl.txt
rm -rf processed_data/ catalog.txt
echo

# Remove downloaded file
rm -rf submitted_data
echo "Tests successful"
