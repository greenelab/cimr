#!/bin/bash

set -e -x

function check_hash {

    data_type=$1

    exdir="examples/expected_processed_data/"${data_type}

    echo "cimr test with data_type ${data_type} example data..."

    cimr processor -process -yaml-file examples/submitted/${data_type}.yml -outdir test

    gzip -dc ${exdir}/${data_type}.txt.gz > ${exdir}/${data_type}.tsv
    expected_hash=$(md5sum ${exdir}/${data_type}.tsv | awk '{print $1}')
    gzip -dc test/${data_type}.txt.gz > test/${data_type}.tsv
    processed_hash=$(md5sum test/${data_type}.tsv | awk '{print $1}')

    if [ "$expected_hash" == "$processed_hash" ]; then
        echo "build test passed."
        rm -rf ${exdir}/${data_type}.tsv
    else
        exit
    fi

}

check_hash gwas

check_hash eqtl

rm -rf test

