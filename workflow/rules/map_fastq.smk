## a snakefile that maps paired FASTQ files indexed by a Biokit sample file
## with bowtie2, filter, and sort, and index the BAM files

## input parameters
##   fq1dict: Dict of forward FASTQ files indexed by sample IDs
##   fq2dict: Dict of reverse FASTQ files indexed by sample IDs
##   bowtie2_index: bowtie2 index of the genome

