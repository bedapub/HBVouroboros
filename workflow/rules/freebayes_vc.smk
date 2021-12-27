import snakemake

sample_annotation = set_samp_anno(False)

samples, fq1dict, fq2dict = parse_sample_annotation(sample_annotation)

rule aggregated_var_inf:
    input:
        expand("results/variant-calling/infref/infref_{sample}_cleaned_allelicprimitives.vcf", sample=samples)
    output:
    	"results/variant-calling/infref/infref_aggregated.vcf"
    shell:
    	"cat {input} > {output}"

rule aggregated_var_inpt:
    input:
        expand("results/variant-calling/inpt/inpt_{sample}_cleaned_allelicprimitives.vcf", sample=samples)
    output:
    	"results/variant-calling/inpt/inpt_aggregated.vcf"
    shell:
    	"cat {input} > {output}"


rule aggregated_var_perSamp:
    input:
        expand("results/variant-calling/perSamp/{sample}/{sample}_cleaned.vcf", sample=samples)
    output:
    	"results/variant-calling/perSamp/perSamp_aggregated.vcf"
    shell:
    	"cat {input} > {output}"


rule correct_bam:
    input:
       refDup ="results/{inpt}/{inpt}_strain.fasta",
       sortBam ="results/{inpt}_bam/{inpt}_{sample}.sorted.bam"
    output:
       "results/{inpt}_bam/{inpt}_{sample}.corrected.sorted.bam"
    shell:
       "workflow/rules/correct_bam.sh -f {input.refDup} -s {input.sortBam} -o {output}"

rule correct_bam_perSamp:
    input:
       refDup ="results/perSamp/{sample}/infref_strain.fasta",
       sortBam ="results/perSamp/{sample}/bam/{sample}.sorted.bam"
    output:
       "results/perSamp/{sample}/bam/{sample}.corrected.sorted.bam"
    shell:
       "workflow/rules/correct_bam.sh -f {input.refDup} -s {input.sortBam} -o {output}"    


rule freebayes_var:
    input:
       refDup ="results/{inpt}/{inpt}_strain_dup.fasta",
       sortBam ="results/{inpt}_bam/{inpt}_{sample}.corrected.sorted.bam"
    output:
       "results/variant-calling/{inpt}/{inpt}_{sample}.vcf"
    shell:
       "freebayes -p 1 -K -m 20 -q 30 --haplotype-length 0 -f  {input.refDup}  {input.sortBam} | vcffilter -f 'QUAL > 20' > {output}"


rule freebayes_var_perSamp:
    input:
       refDup ="results/perSamp/{sample}/infref_strain_dup.fasta",
       sortBam ="results/perSamp/{sample}/bam/{sample}.corrected.sorted.bam"
    output:
       "results/variant-calling/perSamp/{sample}/{sample}.vcf"
    shell:
       "freebayes -p 1 -K -m 20 -q 30  --haplotype-length 0 -f {input.refDup}  {input.sortBam} | vcffilter -f 'QUAL > 20' > {output}"



rule clean_vcfFile:
    input:
       "results/variant-calling/{inpt}/{inpt}_{sample}.vcf"
    output:
       "results/variant-calling/{inpt}/{inpt}_{sample}_cleaned.vcf"
    run:
       vcfClean(str(input),str(output))


rule clean_vcfFile_perSamp:
    input:
       "results/variant-calling/perSamp/{sample}/{sample}.vcf"
    output:
       "results/variant-calling/perSamp/{sample}/{sample}_cleaned.vcf"
    run:
       vcfClean(str(input),str(output))


rule make_vcfallelicprimitives:
    input:
       "results/variant-calling/{inpt}/{inpt}_{sample}_cleaned.vcf"
    output:
       "results/variant-calling/{inpt}/{inpt}_{sample}_cleaned_allelicprimitives.vcf"
    run:
       shell("vcfallelicprimitives {input} > {output}" )


rule make_vcfallelicprimitives_perSamp:
    input:
       "results/variant-calling/perSamp/{sample}/{sample}_cleaned.vcf"
    output:
       "results/variant-calling/perSamp/{sample}/{sample}_cleaned_allelicprimitives.vcf"
    run:
shell("vcfallelicprimitives {input} > {output}" )
