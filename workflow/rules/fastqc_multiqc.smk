"""
Generate fastqc and multiqc reports for FASTQ/BAM files
"""

import os
import pandas as pd
import snakemake

## configfile: "config/config_qc.yaml"

# trimmed files
sample_annotation = config['sample_annotation']

# parse sample annotation
samples, fq1dict, fq2dict = parse_sample_annotation(sample_annotation)

fastqc_dir = "results/fastqc/"
multiqc_dir = "results/multiqc/"
bam_dir = "results/infref_bam/"

rule qc:
    input:
        expand("results/fastqc/{sample}.fastqc.done", sample = samples),
        "results/coverage/infref_genome_depth.done",
        "results/multiqc/multiqc_report.html"

rule fastqc:
    input:
        f1 = lambda wildcards: fq1dict[wildcards.sample],
        f2 = lambda wildcards: fq2dict[wildcards.sample],
    output:
        "results/fastqc/{sample}.fastqc.done"
    shell:
        "fastqc -o {fastqc_dir} {input.f1} {input.f2}; touch {output}"

rule qualimap:
    input:
        "results/infref_bam/{sample}.sorted.bam"
    output:
        temp("results/infref_bam/bamqc/{sample}.bamqc.done")
    shell:
        "qualimap bamqc -bam {input} -c -nw 400 -hm 3 ; touch {output}"

rule covplot:
    input:
        "results/coverage/infref_genome_depth.tsv"
    output:
        temp("results/coverage/infref_genome_depth.done")
    conda:
        "../envs/covplot.yaml"
    shell:
        "Rscript workflow/Rplots.R ; touch {output}"

rule multiqc:
    input:
        expand("results/fastqc/{sample}.fastqc.done",sample = samples),
        expand("results/infref_bam/{sample}.sorted.bam.bai",sample = samples),
        "results/coverage/infref_genome_depth.done"
    output:
        "results/multiqc/multiqc_report.html"
    shell:
        "multiqc {fastqc_dir} {bam_dir} results/coverage/ -o {multiqc_dir}"
