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
cd envs; conda env create; cd -
## in case it has been installed, use the command below to update
## conda env update
conda activate HBVouroboros
```

### Run an example

An out-of-box example can be run by starting the `snakemake` pipeline.

```bash
snakemake -j 99 --use-envmodules ## use --use-conda if no R module is present
```

### Run the pipeline with your own data

Modify the `config/config.yaml` file to specify a sample annotation file.

### Run HBVouroboros using unmapped reads from a Biokit output directory

This feature has been disabled now. It may be activated in the future.

### Validating the sensitivity and specificity of HBVouroboros with RNAsim2

We created RNAsim2, a RNA-seq simulator to validate the sensitivity and
specificity of HBVouroboros. See [RNAsim2/README.md](RNAsim2/README.md) for
details.


## Known issues and solutions

### What to do if conda environment initialization takes too long?

Above we use the default conda solver. If you suffer from slow speed of conda,
consider using [mamba](https://github.com/mamba-org/mamba), which is a drop-in
replacement of conda.

If you met more issues, please raise them using the Issues function of GitHub.
