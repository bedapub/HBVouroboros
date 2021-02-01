import snakemake

if config['doSim'] == True:
	sample_annotation = config['sample_annotation_sm']
else:
	sample_annotation = config['sample_annotation']

samples, fq1dict, fq2dict = parse_sample_annotation(sample_annotation)

rule aggregated_var_inf:
    input:
        expand("results/variant-calling/infref/infref_{sample}_cleaned.vcf", sample=samples)
    output:
    	"results/variant-calling/infref/infref_aggregated.vcf"
    shell:
    	"cat {input} > {output}"

rule aggregated_var_inpt:
    input:
        expand("results/variant-calling/inpt/inpt_{sample}_cleaned.vcf", sample=samples)
    output:
    	"results/variant-calling/inpt/inpt_aggregated.vcf"
    shell:
    	"cat {input} > {output}"

rule freebayes_var:
    input:
       refDup ="results/{inpt}/{inpt}_strain_dup.fasta",
       sortBam ="results/{inpt}_bam/{inpt}_{sample}.sorted.bam"
    output:
       "results/variant-calling/{inpt}/{inpt}_{sample}.vcf"
    shell:
       "freebayes -p 1 -K -m 20 -q 30 -f {input.refDup}  {input.sortBam} > {output}"

rule clean_vcfFile:
    input:
      "results/variant-calling/{inpt}/{inpt}_{sample}.vcf"
    output:
       "results/variant-calling/{inpt}/{inpt}_{sample}_cleaned.vcf"
    run:
       vcfClean(str(input),str(output))


