import snakemake


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
        "workflow/rules/bash/correct_bam.sh -f {input.refDup} -s {input.sortBam} -o {output}"


rule correct_bam_perSamp:
    input:
        refDup ="results/perSamp/{sample}/infref_strain.fasta",
        sortBam ="results/perSamp_bam/{sample}.sorted.bam"
    output:
        "results/perSamp_bam/{sample}.corrected.sorted.bam"
    shell:
        "workflow/rules/bash/correct_bam.sh -f {input.refDup} -s {input.sortBam} -o {output}"    


rule varscan:
    input:
        refDup ="results/{inpt}/{inpt}_strain_dup.fasta",
        sortBam ="results/{inpt}_bam/{inpt}_{sample}.corrected.sorted.bam"
    output:
        pileupFile=temp("results/variant-calling/{inpt}/{inpt}_{sample}_mpileup.tsv"),
        vcfFile="results/variant-calling/{inpt}/{inpt}_{sample}_varscan.vcf"
    shell:
        """
        samtools mpileup -f {input.refDup}  {input.sortBam} >  {output.pileupFile}
        varscan mpileup2snp {output.pileupFile} --min-var-freq 0.01 --output-vcf >  {output.vcfFile}
        """


rule varscan_perSamp:
    input:
        refDup ="results/perSamp/{sample}/infref_strain_dup.fasta",
        sortBam ="results/perSamp_bam/{sample}.corrected.sorted.bam"
    output:
        pileupFile=temp("results/variant-calling/perSamp/{sample}/{sample}_mpileup.tsv"),
        vcfFile="results/variant-calling/perSamp/{sample}/{sample}_varscan.vcf"
    shell:
        """
        samtools mpileup -f {input.refDup}  {input.sortBam} >  {output.pileupFile}
        varscan mpileup2snp {output.pileupFile} --min-var-freq 0.01 --output-vcf >  {output.vcfFile}
        """


rule clean_vcfFile:
    input:
        vcf="results/variant-calling/{inpt}/{inpt}_{sample}_varscan.vcf",
        fasta="results/infref/infref_strain_dup.fasta"
    output:
        "results/variant-calling/{inpt}/{inpt}_{sample}_cleaned.vcf"
    run:
        vcfClean(str(input.vcf), str(input.fasta), str(output))


rule clean_vcfFile_perSamp:
    input:
        vcf="results/variant-calling/perSamp/{sample}/{sample}_varscan.vcf",
        fasta="results/perSamp/{sample}/infref_strain_dup.fasta"
    output:
        "results/variant-calling/perSamp/{sample}/{sample}_cleaned.vcf"
    run:
        vcfClean(str(input.vcf),str(input.fasta),str(output))


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
