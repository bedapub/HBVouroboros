## a snakefile that maps paired FASTQ files indexed by a Biokit sample file
## with bowtie2, filter, and sort, and index the BAM files

## input parameters
##   fq1dict: Dict of forward FASTQ files indexed by sample IDs
##   fq2dict: Dict of reverse FASTQ files indexed by sample IDs
##   bowtie2_index: bowtie2 index of the genome

rule bowtie2_map:
    input:
        f1 = lambda wildcards: fq1dict[wildcards.sample],
        f2 = lambda wildcards: fq2dict[wildcards.sample],
        bowtie2_index = bowtie2_index
    output:
        temp("results/raw_bam/{sample}.bam")
    log:
        "logs/{sample}_bowtie2.log"
    threads:
        8
    shell:
        "bowtie2 -p {threads} --no-mixed --no-discordant --sensitive \
            -x {input.bowtie2_index} \
            -1 {input.f1} -2 {input.f2} 2>{log} | \
            samtools view -Sb - > {output}"


rule filter_and_sort_bam:
    input: "results/raw_bam/{sample}.bam"
    output: "results/bam/{sample}.sorted.bam"
    log:
        "logs/{sample}_filter_and_sort_bam.log"
    threads:
        8
    shell:
        "samtools view -F4 -h {input} | samtools sort -O bam -@ {threads} - > {output}"

rule index_bam:
    input:
        "results/bam/{sample}.sorted.bam"
    output:
        "results/bam/{sample}.sorted.bam.bai"
    threads: 2
    shell:
        "samtools index {input}"

rule flagstat:
    input:
        bam="results/bam/{sample}.sorted.bam",
        bai="results/bam/{sample}.sorted.bam.bai"
    output:
        "results/stats/{sample}.sorted.bam.flagstat"
    threads:
        2
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output}"

