install:
	pip install .

gv: gv/HBVouroboros.gv
	cd gv; $(MAKE)

## run 'ml load pandoc; ml load texlive' first on HPC
README.pdf: README.md gv
	pandoc --number-sections --shift-heading-level-by=-1 -V geometry:"top=3cm, bottom=3cm, left=2cm, right=2cm" README.md -o README.pdf

host=$(hostname)
ifeq "${host}" ""rkalbhpc006""
REFGENOME=/pstore/data/bi/apps/HBVouroboros
else 
REFGENOME=./HBV_refgenomes
endif

HBV_refgenomes:
	bin/HBVouroboros_build_refgenomes "${REFGENOME}"

show_HBV_refgenomes:
	echo "${REFGENOME}"

testSampleAnno=./testdata/sampleAnnotation

test: 
	bin/HBVouroboros_map_samples.py --outdir testdata-out "${REFGENOME}" "${testSampleAnno}"

testBiokitDir=./testbiokit

test_biokit: 
	bin/HBVouroboros_map_biokit.py "${REFGENOME}" "${testBiokitDir}"

.PHONY: install gv test test_biokit
