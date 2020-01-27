#!/bin/bash

set -e -x

cd ~/cimr/.circleci/
date

# default chunksize: 5 million
cimr processor -process -yaml-file dhu_test.yml

# New chunksize: 20 million
#cimr processor -process -yaml-file dhu_test.yml -chunksize 20000000

date
