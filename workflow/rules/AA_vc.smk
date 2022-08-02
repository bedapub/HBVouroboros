import snakemake
# include: "common.smk"
config: "config/config.yaml"

sample_annotation = config['sample_annotation']
samples, fq1dict, fq2dict = parse_sample_annotation(sample_annotation)

#rule all:
#    input:
#        expand("results/variant-calling-AA/{sample}.sorted.sam", sample = samples),
#        expand("results/variant-calling-AA/{sample}.sam2AAFreq.done", sample = samples),
#        "results/infref/inferred_strain_dup.fasta"

rule bam2sam_vc:
    input:
        "results/{inpt}_bam/{inpt}_{sample}.sorted.bam"
    output:
        "results/variant-calling-AA/{inpt}/{inpt}_{sample}.sorted.sam"
    shell:
        "samtools view -h {input} -o {output}"


rule bam2sam_vc_perSamp:
    input:
        "results/perSamp_bam/{sample}.sorted.bam"
    output:
        "results/variant-calling-AA/perSamp/{sample}/{sample}.sorted.sam"
    shell:
        "samtools view -h {input} -o {output}"


rule sam2aaFreq_vc:
    input:
        sam="results/variant-calling-AA/{inpt}/{inpt}_{sample}.sorted.sam",
        fasta= "results/{inpt}/{inpt}_strain_dup.fasta"
    output:
        temp("results/variant-calling-AA/{inpt}/{inpt}_{sample}.sam2AAFreq.done")
    shell:
        "python3 workflow/AA_vc_code/sam2aaFreq.py {input.sam} {input.fasta} ; touch {output}"

rule sam2aaFreq_vc_perSamp:
    input:
        sam="results/variant-calling-AA/perSamp/{sample}/{sample}.sorted.sam",
        fasta= "results/perSamp/{sample}/infref_strain_dup.fasta"
    output:
        temp("results/variant-calling-AA/perSamp/{sample}/{sample}.sam2AAFreq.done")
    shell:
        "python3 workflow/AA_vc_code/sam2aaFreq.py {input.sam} {input.fasta} ; touch {output}"
