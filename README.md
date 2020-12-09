*HBVouroboros* automates sequencing-based HBV genotyping and expression profiling
===

*HBVouroboros* uses RNA-sequencing reads to infer HBV genotype, quantify HBV
transcript expression, and perform variant calling of HBV genomes.

*HBVouroboros*, distributed under the GPL-3 license, is available at
https://github.com/bedapub/HBVouroboros.

## Installation and usage

### Download the source code

```bash
git clone https://github.com/bedapub/HBVouroboros.git
```

### Setup conda environment

```bash
## setup conda environment
cd env; conda env create; cd -
## in case it has been installed, use the command below to update
## conda env update
conda activate HBVouroboros
```

### Run an example

An out-of-box example can be run by starting the `snakemake` pipeline.

```bash
snakemake
```

### Run the pipeline with your own data

Modify the `config/config.yaml` file to specify a sample annotation file.

### Run HBVouroboros using unmapped reads from a Biokit output directory

This feature has been disabled now. It may be activated in the future.

## RNAsim2: a RNA-seq simulator to validate the sensitivity and specificity of
HBVouroboros

An RNA-seq paired-end read simulator is added to the repository in the folder
"RNAsim2". The program outputs two fastq.tar.gz files corresponding to left and
right hand reads, as well as a sample annotation file. These files are deposited
in "RNAsim2/output". In order to run HBVouroboros using the simualted samples,
the config file has to be modified to point to the corresponding sample
annotation file. 

#### Running the RNA-seq simualtor

The program needs to be run from the "RNAsim2/bin" folder. An example call looks
like:

```bash
python RNAsim.py 'gnl|hbvnuc|GQ924620_FT00000_C-C' '45000' '95' '45' --mutate --mutpos "66 88"
```

The above command will generate 45000 paired-end read, with the fragment size of
95 and read size of 45, using the genome  "GQ924620_FT00000_C-C" as reference.
The reference genome will contain random mutations at positions 66 and 88.

