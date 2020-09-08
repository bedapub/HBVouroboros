from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import scripts.refgenomes as ref
from os.path import join

localrules: download_refgenomes

log_dir = "logs"
HTTP = HTTPRemoteProvider()
blastdb_filenames = ["HBV_allgenomes.fasta."+s for s in ("nhr", "nsq", "nin")]

rule all:
    input:
        "HBV_refgenomes_dup_BOWTIE2",
        "HBV_allgenomes.fasta",
        blastdb_filenames

rule download_refgenomes:
    input:
         HTTP.remote('https://hbvdb.lyon.inserm.fr/data/references/hbvdbr.fas',
		     keep_local=True)
    output:
        "HBV_refgenomes.fasta"
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
        "HBV_allgenomes.fasta"
    log:
        join(log_dir, 'download_allgenomes.log')
    run:
        shell("mv {input} {output}")

rule makeblastdb:
    input:
        "HBV_allgenomes.fasta"
    output:
        blastdb_filenames
    log:
        join(log_dir, 'makeblastdb.log')
    shell:
        "makeblastdb -in {input} -title \"HBVdb genomes\" -dbtype nucl 2> {log}"

rule dup_and_conc:
    input:
        fasta="HBV_refgenomes.fasta"
    output:
        fasta="HBV_refgenomes_dup.fasta"
    log:
        join(log_dir, "HBV_refgenomes_dup.log")
    run:
        ref.dup_and_conc_FASTA(input.fasta, output.fasta)

rule bowtie2_index:
    input:
        fasta = "HBV_refgenomes_dup.fasta"
    output:
        directory("HBV_refgenomes_dup_BOWTIE2")
    threads:
        1
    message:
        "Generating bowtie2 index of duplicated HBV reference genomes"
    log:
        join(log_dir, "bowtie2_index.log")
    shell:
        "mkdir {output}; bowtie2-build --threads {threads} {input.fasta} {output}/{output} 2>&1 > {log}; touch {output}/{output}"
