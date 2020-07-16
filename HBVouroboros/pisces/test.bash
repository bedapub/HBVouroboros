#!/bin/bash

if [[ $# != 1 ]]; then
   echo Run the test with $1 OUTDIR
   exit
fi

ml purge
ml load Anaconda3
ml load Singularity/3.5.0

source activate HBVouroboros
snakemake -d $1 --cores 12
