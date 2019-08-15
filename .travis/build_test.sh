#!/bin/bash

set -e -x

cimr processor -process -yaml-file example/submitted/gwas.yml -outdir test
cimr processor -process -yaml-file example/submitted/eqtl.yml -outdir test

rm -rf test

