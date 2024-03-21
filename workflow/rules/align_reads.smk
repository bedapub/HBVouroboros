bowtie2_index = 'resources/ref/HBV_refgenomes_dup_BOWTIE2'
blast_db = 'resources/ref/HBV_allgenomes.fasta'

trinity_fasta = "results/trinity/Trinity.fasta"
trinity_sorted_fasta = "results/trinity/Trinity.sorted.fasta"
blast_out = "results/blast/blast.out"
infref_strain_FASTA = "results/infref/infref_strain.fasta"
infref_strain_dup_FASTA = "results/infref/infref_strain_dup.fasta"
infref_strain_gb = "results/infref/infref_strain.gb"
infref_strain_gff = "results/infref/infref_strain.gff"
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
        "bowtie2 -p {threads} --no-mixed --no-discordant --sensitive -k 1\
            -x {input.bowtie2_index} \
            -1 \"{input.f1}\" -2 \"{input.f2}\" 2>{log} | \
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

rule check_aggregated_fastq:
    input:
        f1 = "results/aggregated_mapped_reads_1.fq.gz",
        f2 = "results/aggregated_mapped_reads_2.fq.gz"
    output:
        "results/stats/aggregated_mapped_reads_count.txt"
    log:
        out = "logs/check_aggregated_fastq.log",
        err = "logs/check_aggregated_fastq.err"
    message:
        "Checking whether aggregated FASTQ files contain valid HBV reads"
    shell:
        "ncount=`zcat results/aggregated_mapped_reads_1.fq.gz | wc -l`; \
         echo \"$ncount\" > {output}; \
         if [ \"$ncount\" -eq \"0\" ]; then \
              echo 'Error: no HBV reads detected. Please check the input.'; \
              echo 'Error: no HBV reads detected. Please check the input.' 2> {log.err} 1> {log.out}; \
              exit 1; \
         fi"

rule run_trinity:
    input:
        f1 = "results/aggregated_mapped_reads_1.fq.gz",
        f2 = "results/aggregated_mapped_reads_2.fq.gz",
        aggcount = "results/stats/aggregated_mapped_reads_count.txt"
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
    input: tsfasta=trinity_sorted_fasta, blastdb={blast_db}
    output: blast_out
    shell:
        "makeblastdb -in {blast_db} -dbtype nucl; blastn -db {blast_db} -query {input.tsfasta} -outfmt 6 > {output}"


rule get_ref_strain_seq_infref:
    input: blast_out
    output: "results/infref/infref_strain.fasta"
    run:
        acc=get_infref_acc(blast_out)
        write_seq_by_acc(blast_db, acc, infref_strain_FASTA)

rule get_ref_strain_gb_infref:
    input: blast_out
    output: "results/infref/infref_strain.gb"
    run:
        gb_acc = get_infref_gb_acc(blast_out)
        download_gb(gb_acc, output[0])


rule ref_strain_gb2gff_inf:
    input:
       "results/infref/infref_strain.gb"
    output:
       "results/infref/infref_strain.gff"
    run:
       gb2gff(input[0], output[0])


rule ref_dup_infref:
    input: "results/infref/infref_strain.fasta"
    output: "results/infref/infref_strain_dup.fasta"
    run:
        dup_and_conc_FASTA(input[0], output[0])


rule bowtie2_index_ref_dup:
    input: "results/{inputRef}/{inputRef}_strain_dup.fasta",
    output: "results/{inputRef}/{inputRef}_bowtie2_index"
    threads: 1
    message:  "Generating bowtie2 index of duplicated infref (and inputRef if present) reference genome"
    log:  "logs/{inputRef}_bowtie2_index_genome.log"
    shell:
        "touch {output}; bowtie2-build --threads {threads} {input} {output}"


rule bowtie2_index_ref_map:
    input:
        genome = 'results/{inputRef}/{inputRef}_bowtie2_index',
        f1 = lambda wildcards: fq1dict[wildcards.sample],
        f2 = lambda wildcards: fq2dict[wildcards.sample]
    output:
        temp("results/{inputRef}_bam/{inputRef}_{sample}.bam")
    log:
        "logs/{inputRef}_bowtie2_{sample}.log"
    threads:
        2
    shell:
        "bowtie2 -p {threads} --no-mixed --no-discordant --sensitive -k 1\
            -x {input.genome} \
            -1 \"{input.f1}\" -2 \"{input.f2}\" 2>{log} | \
            samtools view -Sb - > {output}"


rule filter_and_sort_ref_bam:
    input:
        "results/{inputRef}_bam/{inputRef}_{sample}.bam"
    output:
        "results/{inputRef}_bam/{inputRef}_{sample}.sorted.bam"
    log:
        "logs/{inputRef}_{sample}_filter_and_sort_bam.log"
    threads:
        2
    shell:
        "samtools view -F4 -h {input} | samtools sort -O bam -@ {threads} - > {output}"


rule index_ref_bam:
    input:
        "results/{inputRef}_bam/{inputRef}_{sample}.sorted.bam"
    output:
        "results/{inputRef}_bam/{inputRef}_{sample}.sorted.bam.bai"
    threads: 2
    shell:
        "samtools index {input}"


rule ref_stat:
    input:
        bam = "results/{inputRef}_bam/{inputRef}_{sample}.sorted.bam",
        bai = "results/{inputRef}_bam/{inputRef}_{sample}.sorted.bam.bai"
    output:
        "results/{inputRef}_bam/{inputRef}_{sample}.sorted.bam.stat"
    threads:
        2
    shell:
        "samtools stat -@ {threads} {input.bam} > {output} "


rule genome_count:
    input:
        expand("results/{{inputRef}}_bam/{{inputRef}}_{sample}.sorted.bam.stat",
               sample=samples)
    output:
        "results/coverage/{inputRef}/{inputRef}_genome_count.tsv"
    shell:
         "echo -e 'ID\tcoverage' > {output}; "
         "grep -H '^SN' {input} | \
          grep '1st fragments' | \
          sed 's/.sorted.bam.stat:SN\s1st fragments:\s/\t/g' | \
          sed 's/results\/wildcards.inputRef_bam\///g' >> {output}"





rule dupconc_depth:
    input:
        bam = expand("results/{{inputRef}}_bam/{{inputRef}}_{sample}.sorted.bam",
                     sample=samples),
        bai = expand("results/{{inputRef}}_bam/{{inputRef}}_{sample}.sorted.bam.bai",
                     sample=samples)
    output:
        temp("results/coverage/{inputRef}/{inputRef}_genome_dupconc_depth.tsv")
    shell:
        "samtools depth -a -H -d 0 {input.bam} -o {output} "


rule depth:
    input:
       "results/coverage/{inputRef}/{inputRef}_genome_dupconc_depth.tsv"
    output:
       "results/coverage/{inputRef}/{inputRef}_genome_depth.tsv"
    run:
       dedup_file(input[0], output[0])
       dedup_file(input[0], output[0])


rule dup_gff:
    input:
       gff = "results/{inputRef}/{inputRef}_strain.gff",
       fasta = "results/{inputRef}/{inputRef}_strain_dup.fasta"
    output:
       gff = "results/{inputRef}/{inputRef}_strain_dup.gff"
    run:
       dup_gff(input.fasta, input.gff, output.gff)


rule individual_coverage:
    input:
        bam = "results/{inputRef}_bam/{inputRef}_{sample}.sorted.bam",
        bai = "results/{inputRef}_bam/{inputRef}_{sample}.sorted.bam.bai",
        gff = "results/{inputRef}/{inputRef}_strain_dup.gff"
    output:
        temp("results/coverage/{inputRef}/{inputRef}_genome_{sample}_feature_coverage.tsv")
    shell:
        "coverageBed -counts -a {input.gff} -b {input.bam} > {output}"

## question: can we combine the two?
rule gene_coverage_infref:
     input:
        expand(
          "results/coverage/infref/infref_genome_{sample}_feature_coverage.tsv",
          sample=samples)
     output:
          "results/coverage/infref/infref_genome_gene_coverage.gct"
     run:
        collect_gene_coverage(input, output[0], feat_type='gene')


rule CDS_coverage_infref:
     input:
          expand("results/coverage/infref/infref_genome_{sample}_feature_coverage.tsv", sample=samples)
     output:
          "results/coverage/infref/infref_genome_CDS_coverage.gct"
     run:
          collect_gene_coverage(input, output[0], feat_type='CDS')




################# This section is only executed when "doPerSamp" in config file is set to True



rule aggregate_bam_perSamp:
    	input:
        	"results/bam/{sample}.sorted.bam"
            
    	output:
        	temp("results/bam/aggregated_mapped_reads_{sample}.bam")
    	threads:
        	2
    	shell:
        	"samtools cat --threads {threads} {input} | \
            	samtools sort -n > {output}"


rule perSamp_fq:
    	input:
        	"results/bam/aggregated_mapped_reads_{sample}.bam"
    	output:
        	f1 = temp("results/{sample}_mapped_reads_1.fq.gz"),
        	f2 = temp("results/{sample}_mapped_reads_2.fq.gz")
    	threads: 2
    	shell:
        	"samtools fastq --threads {threads} -N \
            	-1 {output.f1} -2 {output.f2} {input}"


rule run_trinity_perSamp:
    	input:
        	f1 = "results/{sample}_mapped_reads_1.fq.gz",
        	f2 = "results/{sample}_mapped_reads_2.fq.gz",
    	output:
        	"results/perSamp_trinity/{sample}/trinity/Trinity.fasta"
    	threads: 1
    	shell:
        	"Trinity --seqType fq \
            	--left {input.f1} --right {input.f2} \
            	--CPU {threads} --max_memory 10G --output results/perSamp_trinity/{wildcards.sample}/trinity/"


rule run_blast_perSamp:
    	input: "results/perSamp_trinity/{sample}/trinity/Trinity.sorted.fasta"
    	output: "results/perSamp_blast/{sample}_blast.out"
    	shell:
        	"blastn -db {blast_db} -query {input} -outfmt 6 > {output}"


rule sort_trinity_fasta_perSmap:
    	input: "results/perSamp_trinity/{sample}/trinity/Trinity.fasta"
    	output:"results/perSamp_trinity/{sample}/trinity/Trinity.sorted.fasta"
    	run:
        	sort_FASTA_by_length(input[0], output[0])


rule get_ref_strain_seq_perSamp:
    	input: "results/perSamp_blast/{sample}_blast.out"
    	output: "results/perSamp/{sample}/infref_strain.fasta"
    	run:
        	acc=get_infref_acc(input[0])
        	write_seq_by_acc(blast_db, acc, output[0])


rule get_ref_strain_gb_perSamp:
    	input: "results/perSamp_blast/{sample}_blast.out"
    	output: "results/perSamp/{sample}/infref_strain.gb"
    	run:
        	gb_acc = get_infref_gb_acc(input[0])
        	download_gb(gb_acc, output[0])


rule ref_strain_gb2gff_perSamp:
    	input:
       		"results/perSamp/{sample}/infref_strain.gb"
    	output:
       		"results/perSamp/{sample}/infref_strain.gff"
    	run:
       		gb2gff(input[0], output[0])


rule infref_dup_perSamp:
    	input:  "results/perSamp/{sample}/infref_strain.fasta"
    	output: "results/perSamp/{sample}/infref_strain_dup.fasta"
    	run:
        	dup_and_conc_FASTA(input[0], output[0])


rule bowtie2_index_ref_dup_perSmap:
    	input: "results/perSamp/{sample}/infref_strain_dup.fasta"
    	output: "results/perSamp/{sample}/infref_bowtie2_index"
    	threads: 1
	message: "Generating bowtie2 index of duplicated the infref reference genome"
    	log:  "logs/{sample}/bowtie2_index_genome.log"
    	shell: "touch {output}; bowtie2-build --threads {threads} {input} {output}"


rule infref_bowtie2_map_perSamp:
    	input:
        	genome = "results/perSamp/{sample}/infref_bowtie2_index",
        	f1 = lambda wildcards: fq1dict[wildcards.sample],
        	f2 = lambda wildcards: fq2dict[wildcards.sample]
   	output:
        	temp("results/perSamp_bam/{sample}.bam")
	log:
		"logs/{sample}/bowtie2.log"
	threads:
        	2
	shell:
        	"bowtie2 -p {threads} --no-mixed --no-discordant --sensitive -k 1\
		-x {input.genome} -1 \"{input.f1}\" -2 \"{input.f2}\" 2>{log} | \
            	samtools view -Sb - > {output}"


rule filter_and_sort_ref_bam_perSamp:
    	input:
        	"results/perSamp_bam/{sample}.bam"
    	output:
        	"results/perSamp_bam/{sample}.sorted.bam"
    	log:
        	"logs/{sample}/persamp_filter_and_sort_bam.log"
    	threads:
        	2
   	shell:
        	"samtools view -F4 -h {input} | samtools sort -O bam -@ {threads} - > {output}"


rule index_ref_bam_perSamp:
    	input:
       		"results/perSamp_bam/{sample}.sorted.bam"
    	output:
       		"results/perSamp_bam/{sample}.sorted.bam.bai"
    	threads: 2
    	shell:
        	"samtools index {input}"


rule ref_stat_perSamp:
    	input:
        	bam = "results/perSamp_bam/{sample}.sorted.bam",
       		 bai = "results/perSamp_bam/{sample}.sorted.bam.bai"
    	output:
        	"results/perSamp_bam/{sample}.sorted.bam.stat"
    	threads:
        	2
    	shell:
        	"samtools stat -@ {threads} {input.bam} > {output} "


rule genome_count_perSamp:
    	input:
        	"results/perSamp_bam/{sample}.sorted.bam.stat"
    	output:
        	"results/coverage/perSamp/{sample}_genome_count.tsv"
    	shell:
         	"echo -e 'ID\tcoverage' > {output}; "
         	"grep -H '^SN' {input} | \
          	grep '1st fragments' | \
          	sed 's/.sorted.bam.stat:SN\s1st fragments:\s/\t/g' | \
          	sed 's/results\/perSamp_bam\///g' >> {output}"


rule dupconc_depth_perSamp:
    	input:
        	bam = "results/perSamp_bam/{sample}.sorted.bam",
        	bai = "results/perSamp_bam/{sample}.sorted.bam.bai"
    	output:
       		temp("results/coverage/perSamp/{sample}_genome_dupconc_depth.tsv")
    	shell:
        	"samtools depth -a -H -d 0 {input.bam} -o {output}"


rule depth_perSamp:
    	input:
       		"results/coverage/perSamp/{sample}_genome_dupconc_depth.tsv"
    	output:
      		"results/coverage/perSamp/{sample}_genome_depth.tsv"
    	run:
       		dedup_file(input[0], output[0])
       		dedup_file(input[0], output[0])


rule dup_gff_perSamp:
    	input:
       		gff = "results/perSamp/{sample}/infref_strain.gff",
       		fasta = "results/perSamp/{sample}/infref_strain_dup.fasta"
    	output:
       		gff = "results/perSamp/{sample}/infref_strain_dup.gff"
    	run:
       		dup_gff(input.fasta, input.gff, output.gff)


rule individual_coverage_perSamp:
    	input:
      		bam = "results/perSamp_bam/{sample}.sorted.bam",
     		bai = "results/perSamp_bam/{sample}.sorted.bam.bai",
     		gff = "results/perSamp/{sample}/infref_strain_dup.gff"
    	output:
      		"results/coverage/perSamp/{sample}_genome_feature_coverage.tsv"
   	shell:
        	"coverageBed -counts -a {input.gff} -b {input.bam} > {output}"


rule gene_coverage_perSamp:
     	input:
          	"results/coverage/perSamp/{sample}_genome_feature_coverage.tsv"
     	output:
          	"results/coverage/perSamp/{sample}_genome_gene_coverage.gct"
     	run:
        	collect_gene_coverage(input, output[0], feat_type='gene')


rule CDS_coverage_perSamp:
     	input:
        	"results/coverage/perSamp/{sample}_genome_feature_coverage.tsv"
     	output:
          	"results/coverage/perSamp/{sample}_genome_CDS_coverage.gct"
     	run:
          	collect_gene_coverage(input, output[0], feat_type='CDS')

#################


