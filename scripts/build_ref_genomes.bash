#!/bin/bash

# the script builds indices of HBV reference genomes

## the following line is needed to activate conda within a bash script
## see https://github.com/conda/conda/issues/7980
eval "$(conda shell.bash hook)"
conda activate HBVouroboros
snakemake --use-conda -d ../HBV_refgenomes/ $@
