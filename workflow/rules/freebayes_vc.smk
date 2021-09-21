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
    run:
       shell("""set +u; line=$(awk "NR==1{{print $1}}" {input.refDup} | cut -d"=" -f2 | cut -d" " -f1 ); samtools view -h {input.sortBam} | awk -v FS="\\t" -v Line="$line" -v OFS="\\t" '{{{{if ($0 ~ "^[^@]" && $4>Line) {{$4 =$4-Line;}}}} print $0 }}' | samtools sort | samtools view -b  > {output}; set -u """)


rule correct_bam_perSamp:
    input:
       refDup ="results/perSamp/{sample}/infref_strain.fasta",
       sortBam ="results/perSamp/{sample}/bam/{sample}.sorted.bam"
    output:
       "results/perSamp/{sample}/bam/{sample}.corrected.sorted.bam"
    run:
       shell("""set +u; line=$(awk "NR==1{{print $1}}" {input.refDup} | cut -d"=" -f2 | cut -d" " -f1 ); samtools view -h {input.sortBam} | awk -v Line="$line" -v FS="\\t" -v OFS="\\t" '{{{{if ($0 ~ "^[^@]" && $4>Line) {{$4 =$4-Line;}}}} print $0 }}' | samtools sort | samtools view -b > {output}; set -u """)    


rule freebayes_var:
    input:
       refDup ="results/{inpt}/{inpt}_strain_dup.fasta",
       sortBam ="results/{inpt}_bam/{inpt}_{sample}.corrected.sorted.bam"
    output:
       "results/variant-calling/{inpt}/{inpt}_{sample}.vcf"
    shell:
       "freebayes -p 1 -K -m 20 -q 30 -f {input.refDup}  {input.sortBam} > {output}"


rule freebayes_var_perSamp:
    input:
       refDup ="results/perSamp/{sample}/infref_strain_dup.fasta",
       sortBam ="results/perSamp/{sample}/bam/{sample}.corrected.sorted.bam"
    output:
       "results/variant-calling/perSamp/{sample}/{sample}.vcf"
    shell:
       "freebayes -p 1 -K -m 20 -q 30 -f {input.refDup}  {input.sortBam} > {output}"



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



