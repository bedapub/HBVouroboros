## a snakefile that maps paired FASTQ files indexed by a Biokit sample file 
## with bowtie2, filter, and sort, and index the BAM files

## input parameters
## fq1dict: Dict of forward FASTQ files indexed by sample IDs
## fq2dict: Dict of reverse FASTQ files indexed by sample IDs
## bowtie2_index: bowtie2 index of the genome
## rawbam_dir: String, raw bam directory (temporary)
## bam_dir: String, sorted bam directory
## log_dir: String, log directory
## stats_dir: String, stats directory

from os.path import join

rule bowtie2_map:
    input:
        f1 = lambda wildcards: fq1dict[wildcards.sample],
        f2 = lambda wildcards: fq2dict[wildcards.sample],
        bowtie2_index = bowtie2_index,
    output:
        temp(join(rawbam_dir, "{sample}.bam"))
    log:
        join(log_dir, "{sample}_bowtie2.log")
    threads:
        8 
    shell:
        "bowtie2 -p {threads} --no-mixed --no-discordant --sensitive \
            -x {input.bowtie2_index} \
            -1 {input.f1} -2 {input.f2} 2>{log} | \
            samtools view -Sb - > {output}"
    
rule filter_and_sort_bam:
    input:
        join(rawbam_dir, "{sample}.bam")
    output:
        join(bam_dir, "{sample}.sorted.bam")
    log:
        join(log_dir, "{sample}_filter_and_sort_bam.log")
    threads:
        8
    shell:
        "samtools view -F4 -h {input} | samtools sort -O bam -@ {threads} - > {output}"

rule index_bam:
    input:
        join(bam_dir, "{sample}.sorted.bam")
    output:
        join(bam_dir, "{sample}.sorted.bam.bai")
    threads: 2
    shell:
        "samtools index {input}"

rule flagstat:
    input:
        bam=join(bam_dir, "{sample}.sorted.bam"),
        bai=join(bam_dir, "{sample}.sorted.bam.bai"),
    output:
        join(stats_dir, "{sample}.sorted.bam.flagstat")
    threads:
        2
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output}"

