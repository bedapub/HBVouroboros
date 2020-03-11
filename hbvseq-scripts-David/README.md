hbvseq: a bioinformatics software pipeline to handle HBV RNA-sequencing data
===

Contact: Jitao David Zhang <jitao_david.zhang@roche.com>

## Version 1.0

Description: hbvseq is a software pipeline written in Python and R to handle HBV RNA-sequencing data.

Dependencies: python2 (>=2.7) or python 3; emboss (seqret); bowtie; bowtie2; gsnap; samtools (version 1.2); bedtools; trinity

## Version 2.0

We aim at developing a software suite for next-generation sequencing of HBV sequences, based on what we achieved with Version 1.0, however with following improvements in mind

* An intermediate mapping step is not required. Instead, reads are mapped to all available genomes in HBVdb, which are already annotated
* Biokit-like format sample annotation is accepted
* Snakemake is used to orchestrate the pipeline


