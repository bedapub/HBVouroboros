"""
Generate fastqc and multiqc reports for FASTQ/BAM files
"""

import os
import pandas as pd
import snakemake

fastqc_dir = "results/fastqc/"
bam_dir_infref = "results/infref_bam/"
bam_dir_inputRef = "results/inputRef_bam/"
bam_dir_perSamp = "results/perSamp_bam/"

rule qc:
    input:
        expand("results/fastqc/{sample}.fastqc.done", sample = samples),
        "results/coverage/{inputRef}/{inputRef}_genome_depth.done",
        "results/multiqc/multiqc_report.html"

rule fastqc:
    input:
        f1 = lambda wildcards: fq1dict[wildcards.sample],
        f2 = lambda wildcards: fq2dict[wildcards.sample],
    output:
        "results/fastqc/{sample}.fastqc.done"
    shell:
        "fastqc -o results/fastqc/ \"{input.f1}\" \"{input.f2}\"; touch {output}"

rule qualimap:
    input:
        "results/{inputRef}_bam/{inputRef}_{sample}.sorted.bam"
    output:
        temp("results/{inputRef}_bam/bamqc/{sample}.bamqc.done")
    shell:
        "qualimap bamqc -bam {input} -c -nw 400 -hm 3 ; touch {output}"

rule qualimap_persamp:
    input:
        "results/perSamp_bam/{sample}.sorted.bam"
    output:
        temp("results/perSamp_bam/bamqc/{sample}.bamqc.done")
    shell:
        "qualimap bamqc -bam {input} -c -nw 400 -hm 3 ; touch {output}"


rule covplot:
    input:
        "results/coverage/{inputRef}/{inputRef}_genome_depth.tsv"
    output:
        done=temp("results/coverage/{inputRef}/{inputRef}_genome_depth.done"),
        mean="results/coverage/{inputRef}/{inputRef}_genome_depth_mean.tsv", 
        pngf="results/coverage/{inputRef}/{inputRef}_genome_depth_mqc.png"
    shell:
        "Rscript workflow/Rplots.R {input} {output.mean} {output.pngf}; touch {output.done}"


rule covplot_persamp:
    input:
        "results/coverage/perSamp/{sample}_genome_depth.tsv"
    output:
        done=temp("results/coverage/perSamp/{sample}_genome_depth.done"),
        mean="results/coverage/perSamp/{sample}_genome_depth_mean.tsv", 
        pngf="results/coverage/perSamp/{sample}_genome_depth_mqc.png"
    shell:
        "Rscript workflow/Rplots.R {input} {output.mean} {output.pngf}; touch {output.done}"





rule multiqc_inf:
    input:
        fastqcDone=expand("results/fastqc/{sample}.fastqc.done",sample = samples),
        bamBai=expand("results/infref_bam/infref_{sample}.sorted.bam.bai",sample = samples),
        
    output:
        "results/multiqc/infref/infref_multiqc_report.html"
    shell:
        "multiqc --force {fastqc_dir} {bam_dir_infref} results/coverage/infref  --filename 'infref_multiqc_report.html' -o results/multiqc/infref"



rule multiqc_inf_persamp:
    input:
        fastqc=expand("results/fastqc/{sample}.fastqc.done",sample = samples),
        sorted=expand("results/perSamp_bam/{sample}.sorted.bam.bai",sample = samples),
    output:
        "results/multiqc/perSamp/perSamp_multiqc_report.html"
    shell:
        "multiqc -f {fastqc_dir} {bam_dir_perSamp} results/coverage/perSamp  --filename 'perSamp_multiqc_report.html' -o results/multiqc/perSamp"



rule multiqc_inputRef:
    input:
        expand("results/fastqc/{sample}.fastqc.done",sample = samples),
        expand("results/inputRef_bam/inputRef_{sample}.sorted.bam.bai",sample = samples),
    output:
        "results/multiqc/inputRef/inputRef_multiqc_report.html"
    shell:
        "multiqc --force {fastqc_dir} {bam_dir_inputRef} results/coverage/inputRef  --filename 'inputRef_multiqc_report.html' -o results/multiqc/inputRef"
