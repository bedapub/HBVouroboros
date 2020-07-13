import os
import pandas as pd
import csv
from os.path import join, basename, abspath
import snakemake

configfile: snakemake.workflow.srcdir("config/config.yaml")

sample_annotation = config['sample_annotation']
illumina_clip_file = config['illumina_clip_file']
illumina_clip_opts = config['illumina_clip_opts']
trimmomatic_steps = config['trimmomatic_steps']

# parse sample annotation
annotation = pd.read_table(sample_annotation)
samples = annotation.iloc[:, 0]
fq1s = annotation.iloc[:, 2]
fq2s = annotation.iloc[:, 3]

fq1dict = dict(zip(samples, fq1s))
fq2dict = dict(zip(samples, fq2s))

output_sample_annotation_fn="trimmed_" + basename(sample_annotation)
trimmed_dir = abspath("trimmed")
unpaired_dir = abspath("unpaired")
log_dir = abspath("logs")

rule all:
    input:
        output_sample_annotation_fn,
        expand(join(trimmed_dir, "{sample}_R1.fastq.gz"), sample=samples),
        expand(join(trimmed_dir, "{sample}_R2.fastq.gz"), sample=samples),
        expand(join(log_dir, "{sample}_trimmomatic.log.gz"), sample=samples)

rule trimmomatic:
    input:
        f1 = lambda wildcards: abspath(fq1dict[wildcards.sample]),
        f2 = lambda wildcards: abspath(fq2dict[wildcards.sample])
    output:
        tf1 = join(trimmed_dir, "{sample}_R1.fastq.gz"),
        tf2 = join(trimmed_dir, "{sample}_R2.fastq.gz"),
        uf1 = join(unpaired_dir, "{sample}_R1_unpaired.fastq.gz"),
        uf2 = join(unpaired_dir, "{sample}_R2_unpaired.fastq.gz"),
    log:
        logout = join(log_dir, "{sample}_trimmomatic.log"),
        summout = join(log_dir, "{sample}_trimmomatic_summary.txt")
    threads:
        8
    shell:
        "trimmomatic PE -threads {threads} "
        "-trimlog {log.logout} -summary {log.summout} "
        "{input.f1} {input.f2} "
        "{output.tf1} {output.uf1} {output.tf2} {output.uf2} "
        "ILLUMINACLIP:{illumina_clip_file}{illumina_clip_opts} {trimmomatic_steps}"

rule gzip_log:
    input:
        join(log_dir, "{sample}_trimmomatic.log")
    output:
        join(log_dir, "{sample}_trimmomatic.log.gz")
    shell:
        "gzip {input}"

rule output_sample_annotation:
    output:
        tsv = output_sample_annotation_fn
    run:
        outanno = annotation.copy()
        outanno.iloc[:, 2] = [join(trimmed_dir, x+"_R1.fastq.gz") 
            for x in samples]
        outanno.iloc[:, 3] = [join(trimmed_dir, x+"_R2.fastq.gz") 
            for x in samples]
        outanno.to_csv(output.tsv, index=False, sep='\t', quoting=csv.QUOTE_NONE)
