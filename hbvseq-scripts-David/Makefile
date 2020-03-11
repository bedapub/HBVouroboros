PYTHONUSERBASE=/pstore/apps/`id -g -n`/python/`echo ${EBVERSIONPYTHON} |cut -d "." -f 1,2`

test:
	cd hbvseq; nosetests

install: installPy install-R install-sh

installPy:
	pip install --user --upgrade ./hbvseq/

testPy: test.py
	python test.py

install-R:
	cp -p ./rsrc/* ./bin

install-sh:
	cp -p ./src/trinity_from_bams.sh ./bin
