#!/bin/bash

set -ex

cd ~/cimr/tests/gwas-hdl/
cimr processor -process -yaml-file hdl_cholesterol.yml

# Decompress the output file ("-k" option keeps the original gz file).
# Note that WE CAN NOT COMPARE TWO GZ FILES DIRECTLY, because by default a gz file
# includes the original file's name and timestamp in the header when it is generated.
# See more details at: https://github.com/greenelab/cimr-d/issues/14
gunzip -k processed_data/gwas/HDL_Cholesterol.txt.gz

# Verify de-compressed file's md5sum.
# (To disable the output message, use "--status" option.)
# When this command's return status is not zero, CircleCI will quit immediately.
md5sum --check md5.txt
