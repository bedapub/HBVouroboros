import snakemake


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
