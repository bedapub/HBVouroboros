dryrun:
	snakemake -d ../../test-outputs/test-trim-out --dryrun -p

run: clean
	snakemake -d ../../test-outputs/test-trim-out --cores 8 -j 8 --local-cores 8 -p

clean:
	rm -rf ../../test-outputs/test-trim-out

.PHONY=clean
