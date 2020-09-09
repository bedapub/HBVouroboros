*HBVouroboros* automates sequencing-based HBV genotyping and expression profiling
===

*HBVouroboros* uses RNA-sequencing reads to infer HBV genotype, quantify HBV
transcript expression, and perform variant calling of HBV genomes.

*HBVouroboros* is distributed under the GPL-3 license.

The source code is available at https://github.roche.com/BEDA/HBVouroboros.

## The workflow in a nutshell

See Figure 1 for a graphic representation of the workflow.

![The workflow implemented by *HBVouroboros* in a graph](gv/HBVouroboros.svg){ width=300px }

## Installation and usage

### Download the source code

```bash
git clone https://github.com/bedapub/HBVouroboros.git
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
## define the directory of reference genomes here, change the path according
## to your local environment
HBVouroboros_refdir=./HBV_refgenomes
HBVouroboros_build_refgenomes.py "${HBVouroboros_refdir}"
export HBV_refgenomes=${HBVouroboros_refdir}
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
