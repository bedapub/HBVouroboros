import os
import pandas as pd
from os.path import join
from HBVouroboros import biokit
import snakemake

configfile: snakemake.workflow.srcdir("config/config.yaml")

create_genome_size_cmd=config['create_genome_size_cmd']
pisces_cmd=config['pisces_cmd']
pisces_opts=config['pisces_opts']

sample_annotation = config['sample_annotation']
input_fa = config['ref_genome_fa']
ref_genome_str = '{genus} {species} ({build})'.format(
    genus=config['ref_genome_genus'],
    species=config['ref_genome_species'],
    build=config['ref_genome_build'])
    
samples, fq1dict, fq2dict = biokit.parse_sample_annotation(sample_annotation)

ref_genome_dir = "refgenome"
rawbam_dir = "refgenome_rawbam"
bam_dir = "refgenome_bam"
pisces_dir = "refgenome_pisces"
log_dir = "logs"
stats_dir = "stats"
bowtie2_index = join(ref_genome_dir, "genome_bowtie2_index")

include: '../align_reads/map-fastq-snakefile'

ref_genome_GenomeFaFile = join(ref_genome_dir, 'genome.fa')
ref_genome_GenomeSizeFile = join(ref_genome_dir, 'GenomeSize.xml')

rule all:
    input:
        ref_genome_GenomeSizeFile,
        pisces_dir

rule reference_genome_dir:
    output:
        gf=directory(ref_genome_dir)
    shell:
        "mkdir -p {reference_genome_dir}"

rule reference_genome_fa:
    input:
        input_fa
    output:
        ref_genome_GenomeFaFile,
    shell:
        "cp -pr {input} {output}"
        
rule reference_genome_gs:
    input:
        fa=ref_genome_GenomeFaFile
    output:
        gs=ref_genome_GenomeSizeFile
    log:
        join(log_dir, "create-genome-size.log")
    shell:
        "{create_genome_size_cmd} -g {ref_genome_dir} "
        "-s \"{ref_genome_str}\" -o {ref_genome_dir} > {log}"

rule indexing_reference_genome:
    input:
        fa=ref_genome_GenomeFaFile
    output:
        bowtie2_index
    threads: 2
    message: "Generating bowtie2 index of the reference genome"
    log: join(log_dir, "bowtie2-index.log") 
    shell:
        "touch {output}; bowtie2-build --threads {threads} {input} {output}"
        
rule pisces:
    input:
        bam=expand(join(bam_dir, "{sample}.sorted.bam"), sample=samples),
        bai=expand(join(bam_dir, "{sample}.sorted.bam.bai"), sample=samples),
        gs=ref_genome_GenomeSizeFile
    threads: 20
    log: join(log_dir, "pisces.log")
    output:
        directory(pisces_dir)
    shell:
        "{pisces_cmd} --bam {bam_dir} -g {ref_genome_dir} -o {output} "
        "-t {threads} {pisces_opts} > {log}"
