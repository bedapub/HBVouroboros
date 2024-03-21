rule make_report:
    input:
        vcf=expand("results/variant-calling/{{inputRef}}/{{inputRef}}_{sample}_cleaned_allelicprimitives.vcf", sample=samples),
        coverageFiles = "results/coverage/{inputRef}/{inputRef}_genome_count.tsv",
        depth="results/coverage/{inputRef}/{inputRef}_genome_depth.tsv",
        multiqc= "results/multiqc/{inputRef}/{inputRef}_multiqc_report.html", 
        refGenomes="results/{inputRef}/{inputRef}_strain.fasta"

    output:
    	"results/summary/{inputRef}_summary_report.html"
    shell:
    	"python workflow/rules/summary_report.py {wildcards.inputRef} {output}"


rule make_report_persamp:
    input:
        vcf=expand("results/variant-calling/perSamp/{sample}/{sample}_cleaned_allelicprimitives.vcf", sample=samples),
        coverageFiles = expand("results/coverage/perSamp/{sample}_genome_count.tsv", sample=samples),
        depth=expand("results/coverage/perSamp/{sample}_genome_depth.tsv", sample=samples),
        multiqc= "results/multiqc/infref/infref_multiqc_report.html", 
        refGenomes=expand("results/perSamp/{sample}/infref_strain.fasta", sample=samples)

    output:
    	"results/summary/perSamp_summary_report.html"
    shell:
    	"python workflow/rules/summary_report.py 'persamp' {output}"
