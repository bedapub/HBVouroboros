The HBVouroboros workflow
===

# Downloading and building reference genomes

The reference genome of HBV is fetched from a public database called [HBVdb](https://hbvdb.lyon.inserm.fr/HBVdb/HBVdbIndex). 
Previously, the genome was downloaded for every single execution of HBVouroboros via a snakemake file. 
Now, this has been moved to a python script [rules/python/build_refgenomes.py](https://github.com/bedapub/HBVouroboros/blob/code-refactoring/workflow/rules/python/build_refgenomes.py), which only needs to be executed by the user once the data in HBVdb changed.

```bash
python3 rules/python/build_refgenomes.py
```

# Reference

* [Distribution and Reproducibility of Snakemake
  workflows](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)
