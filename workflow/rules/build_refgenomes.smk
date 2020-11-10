from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import refgenomes as ref
from os.path import join

localrules: download_refgenomes

log_dir = "logs"
HTTP = HTTPRemoteProvider()
blastdb_filenames = ["resources/ref/HBV_allgenomes.fasta."+s for s in ("nhr", "nsq", "nin")]

rule download_refgenomes:
    input:
         HTTP.remote('https://hbvdb.lyon.inserm.fr/data/references/hbvdbr.fas',
		     keep_local=True)
    output:
        "resources/ref/HBV_refgenomes.fasta"
    log:
        join(log_dir, 'download_refgenomes.log')
    run:
        shell("mv {input} {output}")

rule download_allgenomes:
    input:
        HTTP.remote(
            'https://hbvdb.lyon.inserm.fr/data/nucleic/fasta/all_Genomes.fas',
            keep_local=True)
    output:
        "resources/ref/HBV_allgenomes.fasta"
    log:
        join(log_dir, 'download_allgenomes.log')
    run:
        shell("mv {input} {output}")

rule makeblastdb:
    input:
        "resources/ref/HBV_allgenomes.fasta"
    output:
        blastdb_filenames
    log:
        join(log_dir, 'makeblastdb.log')
    shell:
        "makeblastdb -in {input} -title \"HBVdb genomes\" -dbtype nucl 2> {log}"

rule dup_and_conc:
    input:
        fasta="resources/ref/HBV_refgenomes.fasta"
    output:
        fasta="resources/ref/HBV_refgenomes_dup.fasta"
    log:
        join(log_dir, "HBV_refgenomes_dup.log")
    run:
        ref.dup_and_conc_FASTA(input.fasta, output.fasta)

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
        join(log_dir, "bowtie2_index.log")
    shell:
        "mkdir {output}; bowtie2-build --threads {threads} {input.fasta} {output} 2>&1 > {log}; touch {output}"
