#!/bin/bash

FASTQPATH=$1
SAMPLEID=$2
OUTPATH=$3

cd $OUTPATH
echo "FASTQS: ${FASTQPATH}"
cellranger vdj --id ${SAMPLEID} \
		--fastqs ${FASTQPATH} \
        --reference=/cellranger/reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \
        --localcores=32 \
		--localmem=30
