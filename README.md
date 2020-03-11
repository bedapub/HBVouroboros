HBVouroboros: A snakemake-based workflow to map RNA-sequencing reads to circular  HBV cccDNA
===
Jitao David Zhang, 11.03.2020

# Background

The ouroboros (or uroboros) is an ancient symbol depicting a serpent or dragon eating its own tail. 

*HBVouroboros* implements a snakemake-based workflow to map Illumina RNA-sequencing reads to HBV cccDNA. *HBVouroboros* offers following functionalities:

* Mapping reads to referece genomes of of eight major HBV genotypes (A-H). The mapping procedure takes care of reads that span the junctions between the two ends of the linear form of the cccDNA.
* Genotype calling
* Base-level and gene-level quantification of read counts
* Reporting

# Usage

# Methods

## HBV reference genomes

We downloaded HBV reference genomes of major genotypes (A-H) from The Hepatitis B Virus Database ([HBVdb](https://hbvdb.lyon.inserm.fr/HBVdb/HBVdbDataset?seqtype=0), release 50, March 2020). HBVdb provides two reference genomes for each of the eight major genotypes, making 16 reference genomes in total. Important information about these reference genomes, including accession number, reference, *etc.*, can be found [in the website of the HBVdb database](https://hbvdb.lyon.inserm.fr/HBVdb/HBVdbNomenclature?nomenclature=table).

## Making aligners work with the circular structure of HBV cccDNA

Most popular sequence aligners that we tested, including bowtie2 and STAR, do not provide specific support for circular genomes. While it is possible to map RNA-sequencing reads derived from HBV cccDNAs to the linear genomes deposited in HBVdb and other public databases, reads that span both ends of the linear genome will not be correctly mapped, causing an inaccurate estimate of HBV gene expression.

We took an approach inspired by the annotation procedure of HBVdb to circumvent this problem[^fn1]. We duplicated the sequence of each reference genome, concatenated the original and the duplicate genome, and used the concatenated genomes (twice the size of the original genome) as templates of read mapping.

[^fn1]: https://academic.oup.com/nar/article/41/D1/D566/1051781

## The mapping procedure

Given that concatenated HBV genomes are used as templates, virtually any sequence alignment software can be used for read mapping. Given that HBV genes are transcribed from cccDNA without splicing, RNA sequencing reads can be aligned in an ungapped fashion. *HBVouroboros* offers the option to choose from two popular aligners: *bowtie* and *STAR*.

## Genotype calling

## Feature count summary

