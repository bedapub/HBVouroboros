""" A snakemake-based workflow to generate fastqc and multiqc reports for FASTQ/BAM files """

import os
import pandas as pd
import snakemake

configfile: snakemake.workflow.srcdir("config/config.yaml")

sample_annotation = config['sample_annotation']

# parse sample annotation
annotation = pd.read_table(sample_annotation)
samples = annotation.iloc[:, 0]
fq1s = annotation.iloc[:, 2]
fq2s = annotation.iloc[:, 3]
fq1dict = dict(zip(samples, fq1s))
fq2dict = dict(zip(samples, fq2s))

fastqc_dir = "results/fastqc/"
multiqc_dir = "results/multiqc/"
bam_dir = "results/bam/"

rule fastqc:
    input:
        f1 = lambda wildcards: fq1dict[wildcards.sample],
        f2 = lambda wildcards: fq2dict[wildcards.sample],
    output:
        os.path.join(fastqc_dir, "{sample}.fastqc.done")
    shell:
        "fastqc -o {fastqc_dir} {input.f1} {input.f2}; touch {output}"

rule multiqc:
    input:
        expand(os.path.join(fastqc_dir, "{sample}.fastqc.done"), sample=samples),
        expand(os.path.join(bam_dir, "{sample}.sorted.bam.bai"), sample=samples)
    output:
        os.path.join(multiqc_dir, "multiqc_report.html")
    shell:
        "multiqc . -o {multiqc_dir}"
