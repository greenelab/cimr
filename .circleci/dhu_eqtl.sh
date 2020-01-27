#!/bin/bash

set -e -x

cd ~/.circleci
date
cimr processor -process -yaml-file dhu_test.yml -chunksize 50000000
date
