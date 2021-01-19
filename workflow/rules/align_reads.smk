import snakemake
import os

include: "common.smk"
configfile: "config/config_qc.yaml"

# trimmed files
sample_annotation = config['sample_annotation']

bowtie2_index = 'resources/ref/HBV_refgenomes_dup_BOWTIE2'
blast_db = 'resources/ref/HBV_allgenomes.fasta'
blastdb_filenames = ["resources/ref/HBV_allgenomes.fasta."+s for s in ("nhr", "nsq", "nin")]

# parse sample annotation
samples, fq1dict, fq2dict = parse_sample_annotation(sample_annotation)

## trinity_outdir="trinity_out_dir"
trinity_fasta = "results/trinity/Trinity.fasta"
trinity_sorted_fasta = "results/trinity/Trinity.sorted.fasta"
blast_out = "results/blast/blast.out"
inferred_strain_FASTA = "results/infref/inferred_strain.fasta"
inferred_strain_dup_FASTA = "results/infref/inferred_strain_dup.fasta"
inferred_strain_gb = "results/infref/inferred_strain.gb"
inferred_strain_gff = "results/infref/inferred_strain.gff"
inferred_strain_dup_gff = "results/infref/inferred_strain_dup.gff"
infref_bowtie2_index = "results/infref/infref_bowtie2_index"

rule all:
    input:
        expand("results/bam/{sample}.bam", sample = samples),
        expand("results/bam/{sample}.sorted.bam",sample = samples),
        expand("results/bam/{sample}.sorted.bam.bai",sample = samples),
        "results/bam/aggregated_mapped_reads.bam",
        "results/aggregated_mapped_reads_1.fq.gz",
        "results/aggregated_mapped_reads_2.fq.gz",
        trinity_fasta,
        trinity_sorted_fasta,
        blast_out,
        inferred_strain_FASTA,
        inferred_strain_gb,
        inferred_strain_gff,
        inferred_strain_dup_FASTA,
        infref_bowtie2_index,
        expand("results/infref_bam/{sample}.bam",sample = samples),
        expand("results/infref_bam/{sample}.temp.bam",sample = samples),
        expand("results/infref_bam/{sample}.temp.bam.bai",sample = samples),
        expand("results/infref_bam/{sample}.nofilter.bam.stat",sample = samples),
        expand("results/infref_bam/{sample}.sorted.bam",sample = samples),
        expand("results/infref_bam/{sample}.sorted.bam.bai",sample = samples),
        expand("results/infref_bam/{sample}.sorted.bam.stat",sample = samples),
        expand("results/coverage/infref_genome_{sample}_feature_coverage.tsv", sample = samples),
        "results/coverage/infref_genome_count.tsv",
        "results/coverage/infref_genome_depth.tsv",
        "results/coverage/infref_genome_gene_coverage.gct",
        "results/coverage/infref_genome_CDS_coverage.gct"

rule bowtie2_map:
    input:
        f1 = lambda wildcards: fq1dict[wildcards.sample],
        f2 = lambda wildcards: fq2dict[wildcards.sample],
        bowtie2_index = bowtie2_index
    output:
        temp("results/bam/{sample}.bam")
    log:
        "logs/{sample}_bowtie2.log"
    threads:
        8
    shell:
        "bowtie2 -p {threads} --no-mixed --no-discordant --sensitive -x {input.bowtie2_index} -1 {input.f1} -2 {input.f2} 2>{log} | samtools view -Sb - > {output}"
        # --un-conc-gz {unmapped_dir}  ## Write paired-end reads that fail to align concordantly to fastq files

rule filter_and_sort_bam:
    input: "results/bam/{sample}.bam"
    output:
        "results/bam/{sample}.sorted.bam",
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
        "samtools flagstat {input.bam} > {output}"

rule aggregate_flagstat:
    input: expand("results/stats/{sample}.sorted.bam.flagstat", sample=samples)
    output: "results/stats/samples.mapping.flagstat"
    shell:
    	"cat {input} > {output}"

rule aggregate_bam:
    input:
        expand("results/bam/{sample}.sorted.bam",
            sample=samples)
    output:
        temp("results/bam/aggregated_mapped_reads.bam")
    threads:
        2
    shell:
        "samtools cat --threads {threads} {input} | \
            samtools sort -n > {output}"

rule aggregate_fq:
    input:
        "results/bam/aggregated_mapped_reads.bam"
    output:
        f1 = temp("results/aggregated_mapped_reads_1.fq.gz"),
        f2 = temp("results/aggregated_mapped_reads_2.fq.gz")
    threads: 2
    shell:
        "samtools fastq --threads {threads} -N \
             -1 {output.f1} -2 {output.f2} {input}"

rule run_trinity:
    input:
        f1 = "results/aggregated_mapped_reads_1.fq.gz",
        f2 = "results/aggregated_mapped_reads_2.fq.gz"
    output:
        trinity_fasta
    threads: 1
    shell:
        "Trinity --seqType fq \
            --left {input.f1} --right {input.f2} \
            --CPU {threads} --max_memory 10G --output results/trinity"

rule sort_trinity_fasta:
    input: trinity_fasta
    output: trinity_sorted_fasta
    run:
        sort_FASTA_by_length(input[0], output[0])

rule run_blast:
    input:
        trinity=trinity_sorted_fasta,
        blasdbfiles=blastdb_filenames
    output: blast_out
    shell:
        "blastn -db {blast_db} -query {input.trinity} -outfmt 6 > {output}"

rule get_ref_strain_seq:
    input: blast_out
    output: inferred_strain_FASTA
    run:
        acc=get_infref_acc(blast_out)
        write_seq_by_acc(blast_db, acc, inferred_strain_FASTA)

rule get_ref_strain_gb:
    input: blast_out
    output: inferred_strain_gb
    run:
        gb_acc = get_infref_gb_acc(blast_out)
        download_gb(gb_acc, output[0])

rule ref_strain_gb2gff:
    input: inferred_strain_gb
    output: inferred_strain_gff
    run:
        gb2gff(input[0], output[0])

rule infref_dup:
    input: inferred_strain_FASTA
    output: inferred_strain_dup_FASTA
    run:
        dup_and_conc_FASTA(input[0], output[0])

rule bowtie2_index_infref_dup:
    input: inferred_strain_dup_FASTA
    output: infref_bowtie2_index
    threads: 1
    message:  "Generating bowtie2 index of duplicated the inferred reference genome"
    log:  "logs/bowtie2_index_infref_genome.log"
    shell:
        "touch {output}; bowtie2-build --threads {threads} {input} {output}"


rule infref_bowtie2_map:
    input:
        genome = infref_bowtie2_index,
        f1 = lambda wildcards: fq1dict[wildcards.sample],
        f2 = lambda wildcards: fq2dict[wildcards.sample]
    output:
        temp("results/infref_bam/{sample}.bam")
    log:
        "logs/{sample}_infref_bowtie2.log"
    threads:
        2
    shell:
        "bowtie2 -p {threads} --no-mixed --no-discordant --sensitive \
            -x {infref_bowtie2_index} \
            -1 {input.f1} -2 {input.f2} 2>{log} | \
            samtools view -Sb - > {output}"

rule sort_infref_bam:
    input:
        "results/infref_bam/{sample}.bam"
    output:
        temp("results/infref_bam/{sample}.temp.bam")
    log:
        "logs/{sample}_infref_temp_bam.log"
    threads:
        2
    shell:
        "samtools sort -O bam -@ {threads} {input} > {output}"

rule index_infref_bam_nofilter:
    input:
        "results/infref_bam/{sample}.temp.bam"
    output:
        temp("results/infref_bam/{sample}.temp.bam.bai")
    threads: 2
    shell:
        "samtools index {input}"

rule infref_stat_nofilter:
    input:
        bam = "results/infref_bam/{sample}.temp.bam",
        bai = "results/infref_bam/{sample}.temp.bam.bai"
    output:
        "results/infref_bam/{sample}.nofilter.bam.stat"
    threads:
        2
    shell:
        "samtools stat -@ {threads} {input.bam} > {output}"        
        
rule filter_and_sort_infref_bam:
    input:
        "results/infref_bam/{sample}.bam"
    output:
        "results/infref_bam/{sample}.sorted.bam"
    log:
        "logs/{sample}_infref_filter_and_sort_bam.log"
    threads:
        2
    shell:
        "samtools view -F4 -h {input} | samtools sort -O bam -@ {threads} - > {output}"

rule index_infref_bam:
    input:
        "results/infref_bam/{sample}.sorted.bam"
    output:
        "results/infref_bam/{sample}.sorted.bam.bai"
    threads: 2
    shell:
        "samtools index {input}"

rule infref_stat:
    input:
        bam = "results/infref_bam/{sample}.sorted.bam",
        bai = "results/infref_bam/{sample}.sorted.bam.bai"
    output:
        "results/infref_bam/{sample}.sorted.bam.stat"
    threads:
        2
    shell:
        "samtools stat -@ {threads} {input.bam} > {output}"

rule genome_count:
    input:
        expand("results/infref_bam/{sample}.sorted.bam.stat",
               sample=samples)
    output:
        "results/coverage/infref_genome_count.tsv"
    shell:
         "echo -e 'ID\tcoverage' > {output}; "
         "grep -H '^SN' {input} | \
          grep '1st fragments' | \
          sed 's/.sorted.bam.stat:SN\s1st fragments:\s/\t/g' | \
          sed 's/results\/infref_bam\///g' >> {output}"

rule dupconc_depth:
    input:
        bam = expand("results/infref_bam/{sample}.sorted.bam",
                     sample=samples),
        bai = expand("results/infref_bam/{sample}.sorted.bam.bai",
                     sample=samples)
    output:
       temp("results/coverage/infref_genome_dupconc_depth.tsv")
    shell:
        "samtools depth -a -H -d 0 {input.bam} -o {output}"

rule depth:
    input:
       "results/coverage/infref_genome_dupconc_depth.tsv"
    output:
       "results/coverage/infref_genome_depth.tsv"
    run:
       dedup_file(input[0], output[0])

rule dup_gff:
    input:
       gff = inferred_strain_gff,
       fasta = inferred_strain_dup_FASTA
    output:
       gff = inferred_strain_dup_gff
    run:
       dup_gff(input.fasta, input.gff, output.gff)

rule individual_coverage:
    input:
        bam = "results/infref_bam/{sample}.sorted.bam",
        bai = "results/infref_bam/{sample}.sorted.bam.bai",
        gff = inferred_strain_dup_gff
    output:
        temp("results/coverage/infref_genome_{sample}_feature_coverage.tsv")
    shell:
        "coverageBed -counts -a {input.gff} -b {input.bam} > {output}"

rule gene_coverage:
     input:
        expand(
          "results/coverage/infref_genome_{sample}_feature_coverage.tsv",
          sample=samples)
     output:
          "results/coverage/infref_genome_gene_coverage.gct"
     run:
        collect_gene_coverage(input, output[0], feat_type='gene')

## question: can we combine the two?
rule CDS_coverage:
     input:
        expand(
          "results/coverage/infref_genome_{sample}_feature_coverage.tsv",
          sample=samples)
     output:
          "results/coverage/infref_genome_CDS_coverage.gct"
     run:
        collect_gene_coverage(input, output[0], feat_type='CDS')

