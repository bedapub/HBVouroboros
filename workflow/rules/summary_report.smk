import snakemake

config: "config/config.yaml"
sample_annotation = config['sample_annotation']

samples, fq1dict, fq2dict = parse_sample_annotation(sample_annotation)


rule make_report:
    input:
        vcf=expand("results/variant-calling/{{inpt}}/{{inpt}}_{sample}_cleaned_allelicprimitives.vcf", sample=samples),
        coverageFiles = "results/coverage/{inpt}/{inpt}_genome_count.tsv",
        depth="results/coverage/{inpt}/{inpt}_genome_depth.tsv",
        multiqc= "results/multiqc/{inpt}/{inpt}_multiqc_report.html", 
        refGenomes="results/{inpt}/{inpt}_strain.fasta"
    output:
    	"results/summary/{inpt}_summary_report.html"
    shell:
    	"python workflow/rules/summary_report.py {wildcards.inpt} {output}"


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
