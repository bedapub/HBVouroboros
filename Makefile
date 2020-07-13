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
testOutDir=./test-outputs
testMapSampleOutDir=${testOutDir}/test-map_samples-out

test: test-map_samples test-biokit test-trimmomatic
test-cluster: test-map_samples-cluster test-biokit-cluster test-trimmomatic-cluster

clean-map_samples:
	rm -rf "${testMapSampleOutDir}"

test-map_samples: clean-map_samples ${HBV_refgenomes} ${testSampleAnno}
	bin/HBVouroboros_map_samples.py --outdir "${testMapSampleOutDir}" "${HBV_refgenomes}" "${testSampleAnno}"  --local

test-map_samples-cluster: clean-map_samples ${HBV_refgenomes} ${testSampleAnno}
	bin/HBVouroboros_map_samples.py --outdir "${testMapSampleOutDir}" "${HBV_refgenomes}" "${testSampleAnno}" 

clean-biokit:
	rm -rf ${testOutDir}/biokit-out

test-biokit: clean-biokit ${HBV_refgenomes} ${testBiokitDir}
	bin/HBVouroboros_map_biokit.py --outdir ${testOutDir}/biokit-out "${HBV_refgenomes}" "${testBiokitDir}" --local

test-biokit-cluster: clean-biokit ${HBV_refgenomes} ${testBiokitDir}
	bin/HBVouroboros_map_biokit.py --outdir ${testOutDir}/biokit-out "${HBV_refgenomes}" "${testBiokitDir}"

clean-trimmomatic:
	rm -rf ${testOutDir}/test-trimmomatic-out

test-trimmomatic: clean-trimmomatic ${testSampleAnno}
	bin/HBVouroboros_trimmomatic.py --outdir ${testOutDir}/test-trimmomatic-out "${testSampleAnno}" --local

test-trimmomatic-cluster: clean-trimmomatic ${testSampleAnno}
	bin/HBVouroboros_trimmomatic.py --outdir ${testOutDir}/test-trimmomatic-out "${testSampleAnno}"

.PHONY: install gv clean-trimmomatic test-map_samples test-biokit test-trimmomatic test-map_samples-cluster test-biokit-cluster test-trimmomatic-cluster clean-map_samples clean-biokit clen-trimmomatic
