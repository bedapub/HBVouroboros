import snakemake

sample_annotation = config['sample_annotation']
samples, fq1dict, fq2dict = parse_sample_annotation(sample_annotation)

rule all:
    input:
        expand("results/Var-calling-results/var-{sample}.vcf",sample=samples)


rule freebayes_var:
     input:
        refDup ="resources/ref/HBV_refgenomes_dup.fasta",
       sortBam = expand("results/bam/{sample}.sorted.bam", sample=samples)
     output:
          "results/Var-calling-results/var-{sample}.vcf"
     shell:
        "freebayes -f {input.refDup}  {input.sortBam} > {output}")
    
