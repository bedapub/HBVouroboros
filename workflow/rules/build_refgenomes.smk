from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from os.path import join

HTTP = HTTPRemoteProvider()
blastdb_filenames = ["resources/ref/HBV_allgenomes.fasta."+s for s in ("nhr", "nsq", "nin")]

## previous rules `download_refgenomes` and `download_allgenomes` are now implemented in bash script download_genomes

rule makeblastdb:
    input:
        "resources/ref/HBV_allgenomes.fasta"
    output:
        blastdb_filenames
    log:
        'logs/ref/makeblastdb.log'
    shell:
        "makeblastdb -in {input} -title \"HBVdb genomes\" -dbtype nucl 2> {log}"

rule dup_and_conc:
    input:
        fasta="resources/ref/HBV_refgenomes.fasta"
    output:
        fasta="resources/ref/HBV_refgenomes_dup.fasta"
    log:
        'logs/ref/HBV_refgenomes_dup.log'
    run:
        dup_and_conc_FASTA(input.fasta, output.fasta)

rule bowtie2_index:
    input:
        fasta = "resources/ref/HBV_refgenomes_dup.fasta"
    output:
        directory("resources/ref/HBV_refgenomes_dup_BOWTIE2")
    threads:
        1
    message:
        "Generating bowtie2 index of duplicated HBV reference genomes"
    log:
        'logs/ref/bowtie2_index.log'
    shell:
        "mkdir {output}; bowtie2-build --threads {threads} {input.fasta} {output} 2>&1 > {log}; touch {output}"
