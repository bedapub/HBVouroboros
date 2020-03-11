Use snakemake to build STAR index of duplicated-concatenated HBV reference genomes
===

## Create conda environment

```bash
conda env create
```

## Activate conda environment

```bash
conda activate HBVouroboros
```

## Build index

```bash
snakemake --use-conda
```
