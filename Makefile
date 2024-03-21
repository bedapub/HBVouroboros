dryrun: workflow/Snakefile
	snakemake -p --dry-run

run: workflow/Snakefile
	snakemake -p --cores 1

clean:
	rm -rf results/*

clean-all-genomes: clean
	rm -rf resources/ref/HBV_allgenomes*

gv: gv/HBVouroboros.gv
	cd gv; $(MAKE)

test:
	pytest -s

## run 'ml load pandoc; ml load texlive' first on HPC
README.pdf: README.md gv
	pandoc --number-sections --shift-heading-level-by=-1 -V geometry:"top=3cm, bottom=3cm, left=2cm, right=2cm" README.md -o README.pdf

.PHONY: gv clean test
