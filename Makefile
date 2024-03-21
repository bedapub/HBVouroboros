dryrun: workflow/Snakefile
	snakemake -p --dry-run --use-envmodules

run: workflow/Snakefile
	snakemake -p -j 99 --use-envmodules

test-run: workflow/Snakemake config/test_config.yaml
	snakemake -p --cores 1 --configfile config/test_config.yaml

pytest:
	pytest -s

clean:
	rm -rf results/*

gv: gv/HBVouroboros.gv
	cd gv; $(MAKE)

## run 'ml load pandoc; ml load texlive' first on HPC
README.pdf: README.md gv
	pandoc --number-sections --shift-heading-level-by=-1 -V geometry:"top=3cm, bottom=3cm, left=2cm, right=2cm" README.md -o README.pdf

.PHONY: gv clean test
