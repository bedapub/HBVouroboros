install:
	pip install .

gv: gv/HBVouroboros.gv
	cd gv; $(MAKE)

## run 'ml load pandoc; ml load texlive' first on HPC
README.pdf: README.md gv
	pandoc --number-sections --shift-heading-level-by=-1 -V geometry:"top=3cm, bottom=3cm, left=2cm, right=2cm" README.md -o README.pdf

ifndef REFGENOME
REFGENOME=/pstore/data/bi/apps/HBVouroboros/
endif

HBV_refgenomes:
	bin/HBVouroboros_build_refgenomes "${REFGENOME}"

show_HBV_refgenomes:
	echo "${REFGENOME}"

testSampleAnno=./testdata/sampleAnnotation
testBiokitDir=./testbiokit

test: ${REFGENOME}
	bin/HBVouroboros_map_samples.py --outdir testdata-out "${REFGENOME}" "${testSampleAnno}"

test_biokit: ${REFGENOME}
	bin/HBVouroboros_map_biokit.py "${REFGENOME}" "${testBiokitDir}"

.PHONY: install gv test test_biokit
