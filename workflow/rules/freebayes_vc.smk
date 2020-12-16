import snakemake

sample_annotation = config['sample_annotation']
samples, fq1dict, fq2dict = parse_sample_annotation(sample_annotation)

rule aggregated_var:
    input:
        expand("results/variant-calling/{sample}.vcf", sample=samples)
    output:
    	"results/variant-calling/aggregated.vcf"
    shell:
    	"cat {input} > {output}"

rule freebayes_var:
    input:
       refDup ="resources/ref/HBV_refgenomes_dup.fasta",
       sortBam = expand("results/bam/{sample}.sorted.bam", sample=samples)
    output:
       "results/variant-calling/{sample}.vcf"
    shell:
        "freebayes -p 1 -K -m 20 -q 30 -f {input.refDup}  {input.sortBam} > {output}"

