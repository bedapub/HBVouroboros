Use snakemake to build STAR index of duplicated-concatenated HBV reference genomes
===

## First-time use: create a conda environment

```bash
conda env create
```

## Subsequent use:

* Option #1: run `build_ref_genomes.bash`
* Option #2: run the following commands explicitly
    1. Activate the conda environment: `conda activate HBVouroboros`
    2. Fetch HBV reference genomes and build index: `snakemake --use-conda -d ../HBV_refgenomes/`
```
