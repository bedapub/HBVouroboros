""" A snakemake-based workflow to generate fastqc and multiqc reports for FASTQ/BAM files """

import os
import pandas as pd
import snakemake

include: "common.smk"
configfile: "config/config_qc.yaml"

# trimmed files
sample_annotation = config['sample_annotation']

# parse sample annotation
samples, fq1dict, fq2dict = parse_sample_annotation(sample_annotation)

fastqc_dir = "results/fastqc/"
multiqc_dir = "results/multiqc/"
bam_dir = "results/stats/"

rule all:
    input:
        expand("results/fastqc/{sample}.fastqc.done", sample = samples),
        "results/multiqc/multiqc_report.html"

rule fastqc:
    input:
        f1 = lambda wildcards: fq1dict[wildcards.sample],
        f2 = lambda wildcards: fq2dict[wildcards.sample],
    output:
        "results/fastqc/{sample}.fastqc.done"
    shell:
        "fastqc -o {fastqc_dir} {input.f1} {input.f2}; touch {output}"

rule multiqc:
    input:
        expand("results/fastqc/{sample}.fastqc.done", sample = samples),
        expand("results/bam/{sample}.sorted.bam.bai", sample = samples)
    output:
        "results/multiqc/multiqc_report.html"
    shell:
        "multiqc {fastqc_dir} {bam_dir} -o {multiqc_dir}"
