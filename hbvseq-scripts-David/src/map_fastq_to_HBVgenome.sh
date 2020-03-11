#!/bin/bash

## Map a series of FASTQ files to the circular HBV genome, report coverage and transcript assembly to identify genotype
## zhangj83 

function help {
    echo "Contact: Jitao David Zhang <jitao_david.zhang@roche.com>, Tel 86251"
}

## parse parameters

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -h|--help)
      help
      exit 0
      ;;
    -s|--samplefile)
      SAMPLEFILE="$2"
      shift
      shift
      ;;
    -r|--refHBVgenome)
      REFGENOME="$2"
      shift
      shift
      ;;
    -o|--outdir)
      OUTDIR="$2"
      shift ## past argument
      shift ## past value
      ;;
done

## Workflow
## (1) If no reference genome is set, reads are mapped to a collection of HBV reference genomes, and Trinity is run over the mapped reads in order to assess the most likely strain
## (2) If reference genome is set, reads are mapped to the extended version of its genome. Trinity is nevertheless run in order to check whether the assembled genome has higher similarity with other genomes. In addition, genome coverage is calculated for each base of the genome.
