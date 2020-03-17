install:
	pip install .

gv: gv/HBVouroboros.gv
	cd gv; $(MAKE)

## run 'ml load pandoc; ml load texlive' first on HPC
README.pdf: README.md gv
	pandoc --number-sections --shift-heading-level-by=-1 -V geometry:"top=3cm, bottom=3cm, left=2cm, right=2cm" README.md -o README.pdf

HBV_refgenomes="`pwd`/HBV_refgenomes"

HBV_refgenomes:
	bin/HBVouroboros_build_refgenomes "${HBV_refgenomes}"

testSampleAnno=./testdata/sampleAnnotation

test: HBV_refgenomes
	bin/HBVouroboros_map_samples.py --outdir testdata-out "${HBV_refgenomes}" "${testSampleAnno}"

.PHONY: install gv
