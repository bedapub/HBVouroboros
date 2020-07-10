*HBVouroboros* automates sequencing-based HBV genotyping and expression profiling
===
[Jitao David Zhang](http://orcid.org/0000-0002-3085-0909)<sup>1</sup>

<sup>1</sup>: Roche Pharma Research and Early Development, Pharmaceutical Sciences, Roche Innovation Center Basel, F. Hoffmann-La Roche Ltd, Grenzacherstrasse 124, Basel, 4070, Switzerland

Corresponding to: [jitao_david.zhang@roche.com](mailto:jitao_david.zhang@roche.com)

## Background

The ouroboros (or uroboros) is an ancient symbol depicting a serpent or dragon eating its own tail. In its twisted form, an ouroboros reassembles covalently closed circular DNA (cccDNA) of hepatitis B virus (HBV).

*HBVouroboros* implements a snakemake-based workflow to map Illumina RNA-sequencing reads to HBV cccDNA. *HBVouroboros* offers following functionalities:

* Mapping reads to reference genomes of of eight major HBV genotypes (A-H). The mapping procedure takes care of reads that span the junctions between the two ends of the linear form of the cccDNA. This holds true also for mapping procedures described later.
* *De novo* assembly of the HBV genome.
* Inference of the reference strain from which the reads are likely generated and genotype calling using BLAST and data from HBVdb.
* Base-level, gene-level, and HBV-genome-level quantification of read counts, as well as structural variants with regard to the inferred reference strain.

## Software availability

*HBVouroboros* is distributed under the GPL-3 license. 

The source code is available at https://github.roche.com/BEDA/HBVouroboros.

## Graphic workflow

See Figure 1 for a graphic representation of the workflow.

![The workflow implemented by *HBVouroboros* in a graph](gv/HBVouroboros.svg){ width=300px }

## Installation

### Download the source code

```bash
git clone git@github.roche.com:BEDA/HBVouroboros.git
```

### Setup conda environment

```bash
## setup conda environment
conda env create 
## in case it has been installed, use the command below to update
## conda env update
conda activate HBVouroboros
```
### Install the python package HBVouroboros in the conda environment

```bash
make install ## alternatively, `pip install .`
```

### Build HBV reference genomes (run once at installation)

```bash
## define the directory of reference genomes here
HBVouroboros_refdir=/pstore/data/bi/apps/HBVouroboros
HBVouroboros_build_refgenomes.py "${HBVouroboros_refdir}"
export HBV_refgenomes=/pstore/data/bi/apps/HBVouroboros/
```

Make sure to export the variable `HBV_refgenomes` before running following commands using Makefile. If not set, the directory `./HBV_refgenomes` is searched in the current path.

## Usage

### Run HBVouroboros using a sample annotation file

```bash
sample_annotation_file=./testdata/sampleAnnotation
HBVouroboros_map_samples.py --outdir testdata-HBVouroboros-outdir ${HBVouroboros_refdir} ${sample_annotation_file}
```

Or, equivalently

```bash
make test ## make sure that the environment variable HBV_refgenomes is set
```

### Run HBVouroboros using unmapped reads from a Biokit output directory

*Biokit* is a software pipeline developed at Bioinformatics and Exploratory Data Analysis group of Roche Pharma Research and Early Development. It exports structured output of RNA-seuquencing read mapping. *HBVouroboros* is able to parse the output directory structure of *Biokit* and run over reads that are not mapped to the host genome (human, mouse, *etc.*). This has the advantage that since the reads mapped to other genomes are filtered, the speed of mapping is much faster than starting from the raw FASTA files.

```bash
biokit_output_dir=~/projects/2020-01-HBVcccDNA-RNAseq/cccDNA_destab_202002/biokit_outdir_cccDNA_destab_PHH_202002
HBVouroboros_map_biokit.py ${HBVouroboros_refdir} ${biokit_output_dir}
```

Or, equivalently

```bash
make test_biokit ## make sure that the environment variable HBV_refgenomes is set
```

Note: from biokit version July 10, 2020, user has to specify the option
`--keep-unmapped` to keep unmapped reads so that HBVouroboros functions
normally.

## Methods

### HBV reference genomes

We downloaded HBV reference genomes of major genotypes (A-H) from The Hepatitis B Virus Database ([HBVdb](https://hbvdb.lyon.inserm.fr/HBVdb/HBVdbDataset?seqtype=0), release 50, March 2020). HBVdb provides two reference genomes for each of the eight major genotypes, making 16 reference genomes in total. Important information about these reference genomes, including accession number, reference, *etc.*, can be found [in the website of the HBVdb database](https://hbvdb.lyon.inserm.fr/HBVdb/HBVdbNomenclature?nomenclature=table).

### Making aligners work with the circular structure of HBV cccDNA

Most popular sequence aligners that we tested, including bowtie2 and STAR, do not provide specific support for circular genomes. While it is possible to map RNA-sequencing reads derived from HBV cccDNAs to the linear genomes deposited in HBVdb and other public databases, reads that span both ends of the linear genome will not be correctly mapped, causing an inaccurate estimate of HBV gene expression.

We took an approach inspired by the annotation procedure of HBVdb to circumvent this problem (https://academic.oup.com/nar/article/41/D1/D566/1051781). We duplicated the sequence of each reference genome, concatenated the original and the duplicate genome, and used the concatenated genomes (twice the size of the original genome) as templates of read mapping.



### Aligning reads to HBV genomes

Given that concatenated HBV genomes are used as templates, virtually any sequence alignment software can be used for read mapping. Given that HBV genes are transcribed from cccDNA without splicing, RNA sequencing reads can be aligned in an ungapped fashion. *HBVouroboros* uses *STA* as the default aligner. Users can modify the Snakemake workflow to use other aligners if it is necessary, for instance to use the same aligner for mapping reads to both HBV and to the host genome or transcriptome, for instance the human genome.

### *De novo* assembly

We use the software *Trinity* to perform *de novo* assembly of the HBV cccDNA genome. The largest contig generated by *Trinity* is used to search against all HBV genomes available in *HBVdb*. The top BLAST hit is used as the *inferred reference strain* for downstream analysis.

### Nucleotide-level, gene-level, and HBV-genome level summary statistics as well as structural variants

*HBVouroboros* reports read coverage per nucleotide base and per gene of the inferred genotype, as well as single-nucleotide polymorphisms (SNP) compared with the inferred reference genome. 

## Acknowledgement

I thank Roland Ambs, Raphael Sebastian Müller, and Alex Seitz for sharing their experience with Snakemake. I thank Yaling Zhang, Vivian Wang, Henrik Müller, Miriam Triyatni, Birian Leonard, Souphalone Luangsay, and other colleagues who provided datasets on which HBVouroboros was applied, and for giving feedbacks which allowed HBVouroboros to be continuously improved. I thank Roland Schmucki and Fabian Birzele for their work on the Biokit pipeline.


