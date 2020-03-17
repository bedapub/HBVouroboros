install:
	pip install .

gv:
	cd gv; $(MAKE)

## run 'ml load pandoc' first on HPC
README.pdf: README.md
	pandoc README.md -o README.pdf

.PHONY: install gv
