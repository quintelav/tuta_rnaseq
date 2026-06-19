#!/bin/bash

PROJECT_DIR=~/projects/tuta_rnaseq

mkdir -p $PROJECT_DIR/raw_data

while read SRR
do
    echo "==================================="
    echo "Downloading $SRR"
    echo "==================================="

    prefetch $SRR

    fasterq-dump $SRR \
        --split-files \
        -O $PROJECT_DIR/raw_data \
        -t $PROJECT_DIR/tmp\
        -e 6

    gzip $PROJECT_DIR/raw_data/${SRR}_*.fastq

    echo "$SRR completed"

done < $PROJECT_DIR/metadata/test_samples.txt
