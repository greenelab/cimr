#!/bin/bash

set -e -x

cd ~/cimr/.circleci/
date

# default chunksize: 5 million
#cimr processor -process -yaml-file dhu_test.yml

# New chunksize: 1 million
cimr processor -process -yaml-file dhu_test.yml -chunksize 1000000

date
