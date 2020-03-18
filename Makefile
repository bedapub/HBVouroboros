install:
	pip install .

gv: gv/HBVouroboros.gv
	cd gv; $(MAKE)

## run 'ml load pandoc; ml load texlive' first on HPC
README.pdf: README.md gv
	pandoc --number-sections --shift-heading-level-by=-1 -V geometry:"top=3cm, bottom=3cm, left=2cm, right=2cm" README.md -o README.pdf

ifndef HBV_refgenomes
HBV_refgenomes=./HBV_refgenomes
endif

${HBV_refgenomes}:
	bin/HBVouroboros_build_refgenomes.py "${HBV_refgenomes}"

show_HBV_refgenomes:
	echo "${HBV_refgenomes}"

testSampleAnno=./testdata/sampleAnnotation
testBiokitDir=./testbiokit

test: ${HBV_refgenomes} ${testSampleAnno}
	bin/HBVouroboros_map_samples.py --outdir testdata-out "${HBV_refgenomes}" "${testSampleAnno}"  --local

test-cluster: ${HBV_refgenomes} ${testSampleAnno}
	bin/HBVouroboros_map_samples.py --outdir testdata-out "${HBV_refgenomes}" "${testSampleAnno}" 

test_biokit: ${HBV_refgenomes} ${testBiokitDir}
	bin/HBVouroboros_map_biokit.py "${HBV_refgenomes}" "${testBiokitDir}" --local

test_biokit-cluster: ${HBV_refgenomes} ${testBiokitDir}
	bin/HBVouroboros_map_biokit.py "${HBV_refgenomes}" "${testBiokitDir}"

.PHONY: install gv test test_biokit
