#!/bin/bash

FASTQPATH=$1
SAMPLEID=$2
OUTPATH=$3

cd $OUTPATH
echo "FASTQS: ${FASTQPATH}"
cellranger count --id ${SAMPLEID} \
                --fastqs ${FASTQPATH} \
                --transcriptome=/cellranger/reference/refdata-gex-mm10-2020-A/ \
                --localcores=32 \
                 --localmem=30
