The HBVouroboros workflow
===

# Downloading and building reference genomes

Now we expose the Snakefile directly to users, which downloads and builds
indexes of HBV genomes. This step was done by the wrapper
`bin/HBVouroboros_build_refgenomes.py` before. The implementation now does not
need a wrapper, and the users can control the details of building the genome.

```bash
snakemake --snakefile build_refgenomes.snakefile --cores 1 --directory
../results/testHBVgenome
```

# Reference

* [Distribution and Reproducibility of Snakemake
  workflows](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)
