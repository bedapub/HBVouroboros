import snakemake

if config['doSim'] == True:
	sample_annotation = config['sample_annotation_sm']
else:
	sample_annotation = config['sample_annotation']

bowtie2_index = 'resources/ref/HBV_refgenomes_dup_BOWTIE2'
blast_db = 'resources/ref/HBV_allgenomes.fasta'

# parse sample annotation
samples, fq1dict, fq2dict = parse_sample_annotation(sample_annotation)

## trinity_outdir="trinity_out_dir"
trinity_fasta = "results/trinity/Trinity.fasta"
trinity_sorted_fasta = "results/trinity/Trinity.sorted.fasta"
blast_out = "results/blast/blast.out"
infref_strain_FASTA = "results/infref/infref_strain.fasta"
inpt_strain_FASTA = "results/inpt/inpt_strain.fasta"
infref_strain_gb = "results/infref/infref_strain.gb"
inpt_strain_gb = "results/inpt/inpt_strain.gb"
infref_strain_gff = "results/infref/infref_strain.gff"
inpt_strain_gff = "results/inpt/inpt_strain.gff"
infref_strain_dup_gff = "results/infref/infref_strain_dup.gff"
infref_bowtie2_index = "results/infref/infref_bowtie2_index"

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
        "bowtie2 -p {threads} --no-mixed --no-discordant --sensitive \
            -x {input.bowtie2_index} \
            -1 {input.f1} -2 {input.f2} 2>{log} | \
            samtools view -Sb - > {output}"


rule filter_and_sort_bam:
    input: "results/bam/{sample}.bam"
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
        "results/bam/aggregated_mapped_reads.bam"
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

rule run_trinity_perSamp:
    input:
        f1 = lambda wildcards: fq1dict[wildcards.sample],
        f2 = lambda wildcards: fq2dict[wildcards.sample],
    output:
        "results/perSamp_trinity/{sample}_Trinity.fasta"
    threads: 1
    shell:
        "Trinity --seqType fq \
            --left {input.f1} --right {input.f2} \
            --CPU {threads} --max_memory 10G --output {output}"

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


rule sort_trinity_fasta_perSmap:
    input: "results/perSamp_trinity/{sample}_Trinity.fasta"
    output: "results/perSamp_trinity/{sample}_Trinity.sorted.fasta"
    run:
        sort_FASTA_by_length(input, output)

rule sort_trinity_fasta:
    input: trinity_fasta
    output: trinity_sorted_fasta
    run:
        sort_FASTA_by_length(input[0], output[0])


rule run_blast_perSamp:
    input: "results/perSamp_trinity/{sample}_Trinity.sorted.fasta"
    output: "results/persamp_blast/{sample}_blast.out"
    shell:
        "blastn -db {blast_db} -query {input} -outfmt 6 > {output}"

rule run_blast:
    input: trinity_sorted_fasta
    output: blast_out
    shell:
        "blastn -db {blast_db} -query {input} -outfmt 6 > {output}"

acc_inpt = config['inputRef']
gb_acc_inpt = acc_inpt.split("|")[2].split("_")[0]

rule get_ref_strain_seq_inpt:
    input: blast_out
    output: "results/inpt/inpt_strain.fasta"
    run:
        write_seq_by_acc(blast_db, acc_inpt, inpt_strain_FASTA)

rule get_ref_strain_gb_inpt:
    input: blast_out
    output: "results/inpt/inpt_strain.gb"
    run:
        download_gb(gb_acc_inpt, output[0])

rule get_ref_strain_seq:
    input: blast_out
    output: "results/infref/infref_strain.fasta"
    run:
        acc=get_infref_acc(blast_out)
        write_seq_by_acc(blast_db, acc, infref_strain_FASTA)

rule get_ref_strain_gb:
    input: blast_out
    output: "results/infref/infref_strain.gb"
    run:
        gb_acc = get_infref_gb_acc(blast_out)
        download_gb(gb_acc, output[0])

rule ref_strain_gb2gff:
    input:
       "results/infref/infref_strain.gb", "results/inpt/inpt_strain.gb"
    output:
       "results/infref/infref_strain.gff", "results/inpt/inpt_strain.gff"
    run:
       gb2gff(input[0], output[0])
       gb2gff(input[1], output[1])


rule infref_dup:
    input: expand("results/{inpt}/{inpt}_strain.fasta", inpt=["inpt","infref"])
    output: expand("results/{inpt}/{inpt}_strain_dup.fasta", inpt=["inpt","infref"])
    run:
        dup_and_conc_FASTA(input[0], output[0])
	dup_and_conc_FASTA(input[1], output[1])
	
rule bowtie2_index_infref_dup:
    input: "results/{inpt}/{inpt}_strain_dup.fasta",
    output: "results/{inpt}/{inpt}_bowtie2_index"
    threads: 1
    message:  "Generating bowtie2 index of duplicated the infref reference genome"
    log:  "logs/{inpt}_bowtie2_index_genome.log"
    shell:
        "touch {output}; bowtie2-build --threads {threads} {input} {output}"



rule infref_bowtie2_map:
    input:
        genome = 'results/{inpt}/{inpt}_bowtie2_index',	
        f1 = lambda wildcards: fq1dict[wildcards.sample],
        f2 = lambda wildcards: fq2dict[wildcards.sample]
    output:
        temp("results/{inpt}_bam/{inpt}_{sample}.bam")
    log:
        "logs/{inpt}_bowtie2_{sample}.log"
    threads:
        2
    shell:
        "bowtie2 -p {threads} --no-mixed --no-discordant --sensitive \
            -x {input.genome} \
            -1 {input.f1} -2 {input.f2} 2>{log} | \
            samtools view -Sb - > {output}"

rule filter_and_sort_infref_bam:
    input:
        "results/{inpt}_bam/{inpt}_{sample}.bam"
    output:
        "results/{inpt}_bam/{inpt}_{sample}.sorted.bam"
    log:
        "logs/{inpt}_{sample}_infref_filter_and_sort_bam.log"
    threads:
        2
    shell:
        "samtools view -F4 -h {input} | samtools sort -O bam -@ {threads} - > {output}"

rule index_infref_bam:
    input:
        "results/{inpt}_bam/{inpt}_{sample}.sorted.bam"
    output:
        "results/{inpt}_bam/{inpt}_{sample}.sorted.bam.bai"
    threads: 2
    shell:
        "samtools index {input}"

rule infref_stat:
    input:
        bam = "results/{inpt}_bam/{inpt}_{sample}.sorted.bam",
        bai = "results/{inpt}_bam/{inpt}_{sample}.sorted.bam.bai"
    output:
        "results/{inpt}_bam/{inpt}_{sample}.sorted.bam.stat"
    threads:
        2
    shell:
        "samtools stat -@ {threads} {input.bam} > {output} "

rule genome_count_infref:
    input:
        expand("results/infref_bam/infref_{sample}.sorted.bam.stat",
               sample=samples)
    output:
        "results/coverage/infref_genome_count.tsv"
    shell:
         "echo -e 'ID\tcoverage' > {output}; "
         "grep -H '^SN' {input} | \
          grep '1st fragments' | \
          sed 's/.sorted.bam.stat:SN\s1st fragments:\s/\t/g' | \
          sed 's/results\/infref_bam\///g' >> {output}"

rule genome_count_inpt:
    input:
        expand("results/inpt_bam/inpt_{sample}.sorted.bam.stat",
               sample=samples)
    output:
        "results/coverage/inpt_genome_count.tsv"
    shell:
         "echo -e 'ID\tcoverage' > {output}; "
         "grep -H '^SN' {input} | \
          grep '1st fragments' | \
          sed 's/.sorted.bam.stat:SN\s1st fragments:\s/\t/g' | \
          sed 's/results\/inpt_bam\///g' >> {output}"


rule dupconc_depth_inpt:
    input:
        bam = expand("results/inpt_bam/inpt_{sample}.sorted.bam",
                     sample=samples),
        bai = expand("results/inpt_bam/inpt_{sample}.sorted.bam.bai",
                     sample=samples)
    output:
        temp("results/coverage/inpt_genome_dupconc_depth.tsv")
    shell:
        "samtools depth -a -H -d 0 {input.bam} -o {output} "


rule dupconc_depth_infref:
    input:
        bam = expand("results/infref_bam/infref_{sample}.sorted.bam",
                     sample=samples),
        bai = expand("results/infref_bam/infref_{sample}.sorted.bam.bai",
                     sample=samples)
    output:
       temp("results/coverage/infref_genome_dupconc_depth.tsv")
    shell:
        "samtools depth -a -H -d 0 {input.bam} -o {output}"

rule depth:
    input:
       "results/coverage/{inpt}_genome_dupconc_depth.tsv"
    output:
       "results/coverage/{inpt}_genome_depth.tsv"
    run:
       dedup_file(input[0], output[0])
       dedup_file(input[0], output[0])

rule dup_gff:
    input:
       gff = "results/{inpt}/{inpt}_strain.gff",
       fasta = "results/{inpt}/{inpt}_strain.fasta"
    output:
       gff = "results/{inpt}/{inpt}_strain_dup.gff"
    run:
       dup_gff(input.fasta, input.gff, output.gff)

rule individual_coverage:
    input:
        bam = "results/{inpt}_bam/{inpt}_{sample}.sorted.bam",
        bai = "results/{inpt}_bam/{inpt}_{sample}.sorted.bam.bai",
        gff = "results/{inpt}/{inpt}_strain_dup.gff"
    output:
        temp("results/coverage/{inpt}_genome_{sample}_feature_coverage.tsv")
    shell:
        "coverageBed -counts -a {input.gff} -b {input.bam} > {output}"

rule gene_coverage:
     input:
        input=expand(
          "results/coverage/inpt_genome_{sample}_feature_coverage.tsv",
          sample=samples),
        infref=expand(
          "results/coverage/infref_genome_{sample}_feature_coverage.tsv",
          sample=samples)
     output:
          input="results/coverage/inpt_genome_gene_coverage.gct",
	  infref="results/coverage/infref_genome_gene_coverage.gct"
     run:
        collect_gene_coverage(input.input, output.input, feat_type='gene')
	collect_gene_coverage(input.infref, output.infref, feat_type='gene')

## question: can we combine the two?
rule CDS_coverage:
     input:
        input=expand(
          "results/coverage/inpt_genome_{sample}_feature_coverage.tsv",
          sample=samples),
        inf=expand(
          "results/coverage/infref_genome_{sample}_feature_coverage.tsv",
          sample=samples)
     output:
          input="results/coverage/inpt_genome_CDS_coverage.gct",
	  inf="results/coverage/infref_genome_CDS_coverage.gct"
     run:
        collect_gene_coverage(input.inf, output.inf, feat_type='CDS')
	collect_gene_coverage(input.input, output.input, feat_type='CDS')
