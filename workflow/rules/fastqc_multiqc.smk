"""
Generate fastqc and multiqc reports for FASTQ/BAM files
"""

import os
import pandas as pd
import snakemake

config: "config/config_qc.yaml"

# trimmed files
sample_annotation = config['sample_annotation']

# parse sample annotation
samples, fq1dict, fq2dict = parse_sample_annotation(sample_annotation)

fastqc_dir = "results/fastqc/"
bam_dir_infref = "results/infref_bam/"
bam_dir_inpt = "results/inpt_bam/"

rule qc:
    input:
        expand("results/fastqc/{sample}.fastqc.done", sample = samples),
        "results/coverage/{inpt}/{inpt}_genome_depth.done",
        "results/multiqc/multiqc_report.html"

rule fastqc:
    input:
        f1 = lambda wildcards: fq1dict[wildcards.sample],
        f2 = lambda wildcards: fq2dict[wildcards.sample],
    output:
        "results/fastqc/{sample}.fastqc.done"
    shell:
        "fastqc -o results/fastqc/ {input.f1} {input.f2}; touch {output}"

rule qualimap:
    input:
        "results/{inpt}_bam/{inpt}_{sample}.sorted.bam"
    output:
        temp("results/{inpt}_bam/bamqc/{sample}.bamqc.done")
    shell:
        "qualimap bamqc -bam {input} -c -nw 400 -hm 3 ; touch {output}"

rule covplot:
    input:
        "results/coverage/{inpt}/{inpt}_genome_depth.tsv"
    output:
        done=temp("results/coverage/{inpt}/{inpt}_genome_depth.done"),
        mean="results/coverage/{inpt}/{inpt}_genome_depth_mean.tsv", 
        pngf="results/coverage/{inpt}/{inpt}_genome_depth_mqc.png"
    shell:
        "Rscript workflow/Rplots.R {input} {output.mean} {output.pngf}; touch {output.done}"


rule multiqc_inf:
    input:
        expand("results/fastqc/{sample}.fastqc.done",sample = samples),
        expand("results/infref_bam/infref_{sample}.sorted.bam.bai",sample = samples),
    output:
        "results/multiqc/infref/infref_multiqc_report.html"
    shell:
        "multiqc -f {fastqc_dir} {bam_dir_infref} results/coverage/infref  --filename 'infref_multiqc_report.html' -o results/multiqc/infref"

rule multiqc_inpt:
    input:
        expand("results/fastqc/{sample}.fastqc.done",sample = samples),
        expand("results/inpt_bam/inpt_{sample}.sorted.bam.bai",sample = samples),
    output:
        "results/multiqc/inpt/inpt_multiqc_report.html"
    shell:
        "multiqc -f {fastqc_dir} {bam_dir_inpt} results/coverage/inpt  --filename 'inpt_multiqc_report.html' -o results/multiqc/inpt"
