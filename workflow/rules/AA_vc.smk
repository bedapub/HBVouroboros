import snakemake
# include: "common.smk"
# configfile: "config/config.yaml"

sample_annotation = config['sample_annotation']
samples, fq1dict, fq2dict = parse_sample_annotation(sample_annotation)

#rule all:
#    input:
#        expand("results/variant-calling-AA/{sample}.sorted.sam", sample = samples),
#        expand("results/variant-calling-AA/{sample}.sam2AAFreq.done", sample = samples),
#        "results/infref/inferred_strain_dup.fasta"

rule bam2sam_vc:
    input:
        "results/infref_bam/{sample}.sorted.bam"
    output:
        "results/variant-calling-AA/{sample}.sorted.sam"
    shell:
        "samtools view -h {input} -o {output}"

rule sam2aaFreq_vc:
    input:
        "results/variant-calling-AA/{sample}.sorted.sam"
    output:
        temp("results/variant-calling-AA/{sample}.sam2AAFreq.done")
    shell:
        "python3 workflow/AA_vc_code/sam2aaFreq.py {input}; touch {output}"
