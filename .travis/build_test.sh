#!/bin/bash

set -e -x

cimr processor -process -yaml-file examples/expected_submitted/gwas.yml -outdir test
cimr processor -process -yaml-file examples/expected_submitted/eqtl.yml -outdir test

rm -rf test

