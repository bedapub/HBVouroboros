install:
	pip install .

gv:
	cd gv; $(MAKE)

.PHONY: install gv
