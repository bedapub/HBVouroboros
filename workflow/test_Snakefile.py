import os
import pytest
import argparse
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from subprocess import Popen, PIPE, STDOUT

sys.path.insert(0, 'workflow/rules/')
import common

#include: "rules/common.py"

def test_snakemake_output_files():
    #Run the pipeline
    p= Popen('snakemake --cores 20 --use-conda --latency-wait 900', shell=True, stdout= PIPE, stderr= STDOUT)
    pout= p.stdout.read()
    print(pout.decode('utf-8'))
     
    #Optionally, get stdout and stderr
    stdout, stderr= p.communicate()

    #Check exits code and other expected output            
    #assert 0 == p.returncode

    #Read bam files
    assert os.path.exists('results/bam/simSample1.sorted.bam') == 1
    assert os.path.exists('results/bam/simSample2.sorted.bam') == 1

    #BAM and related files
    assert os.path.exists('results/infref_bam/infref_simSample1.sorted.bam.stat') == 1
    assert os.path.exists('results/infref_bam/infref_simSample1.sorted.bam') == 1
    assert os.path.exists('results/infref_bam/infref_simSample1.corrected.sorted.bam') == 1
    assert os.path.exists('results/infref_bam/infref_simSample1.sorted.bam.bai') == 1
    assert os.path.exists('results/infref_bam/infref_simSample2.sorted.bam.stat') == 1
    assert os.path.exists('results/infref_bam/infref_simSample2.sorted.bam') == 1
    assert os.path.exists('results/infref_bam/infref_simSample2.corrected.sorted.bam') == 1
    assert os.path.exists('results/infref_bam/infref_simSample2.sorted.bam.bai') == 1
    assert os.path.exists('results/inpt_bam/inpt_simSample1.sorted.bam.stat') == 1
    assert os.path.exists('results/inpt_bam/inpt_simSample1.sorted.bam') == 1
    assert os.path.exists('results/inpt_bam/inpt_simSample1.corrected.sorted.bam') == 1
    assert os.path.exists('results/inpt_bam/inpt_simSample1.sorted.bam.bai') == 1
    assert os.path.exists('results/inpt_bam/inpt_simSample2.sorted.bam.stat') == 1
    assert os.path.exists('results/inpt_bam/inpt_simSample2.sorted.bam') == 1
    assert os.path.exists('results/inpt_bam/inpt_simSample2.corrected.sorted.bam') == 1
    assert os.path.exists('results/inpt_bam/inpt_simSample2.sorted.bam.bai') == 1
    assert os.path.exists('results/perSamp/simSample1/bam/simSample1.sorted.bam.stat') == 1
    assert os.path.exists('results/perSamp/simSample1/bam/simSample1.sorted.bam') == 1
    assert os.path.exists('results/perSamp/simSample1/bam/simSample1.corrected.sorted.bam') == 1
    assert os.path.exists('results/perSamp/simSample1/bam/simSample1.sorted.bam.bai') == 1
    assert os.path.exists('results/perSamp/simSample2/bam/simSample2.sorted.bam.stat') == 1
    assert os.path.exists('results/perSamp/simSample2/bam/simSample2.sorted.bam') == 1
    assert os.path.exists('results/perSamp/simSample2/bam/simSample2.corrected.sorted.bam') == 1
    assert os.path.exists('results/perSamp/simSample2/bam/simSample2.sorted.bam.bai') == 1

    #Strain fasta, gff, gb and bowtie index files
    assert os.path.exists('results/infref/infref_strain.gb') == 1
    assert os.path.exists('results/infref/infref_strain_dup.fasta') == 1
    assert os.path.exists('results/infref/infref_strain_dup.gff') == 1
    assert os.path.exists('results/infref/infref_strain.fasta') == 1
    assert os.path.exists('results/infref/infref_bowtie2_index') == 1
    assert os.path.exists('results/infref/infref_strain.gff') == 1
    assert os.path.exists('results/inpt/inpt_strain.gb') == 1
    assert os.path.exists('results/inpt/inpt_strain_dup.fasta') == 1
    assert os.path.exists('results/inpt/inpt_strain_dup.gff') == 1
    assert os.path.exists('results/inpt/inpt_strain.fasta') == 1
    assert os.path.exists('results/inpt/inpt_bowtie2_index') == 1
    assert os.path.exists('results/inpt/inpt_strain.gff') == 1
    assert os.path.exists('results/perSamp/simSample1/infref_strain.gb') == 1
    assert os.path.exists('results/perSamp/simSample1/infref_strain_dup.fasta') == 1
    assert os.path.exists('results/perSamp/simSample1/infref_strain_dup.gff') == 1
    assert os.path.exists('results/perSamp/simSample1/infref_strain.fasta') == 1
    assert os.path.exists('results/perSamp/simSample1/infref_bowtie2_index') == 1
    assert os.path.exists('results/perSamp/simSample1/infref_strain.gff') == 1
    assert os.path.exists('results/perSamp/simSample2/infref_strain.gb') == 1
    assert os.path.exists('results/perSamp/simSample2/infref_strain_dup.fasta') == 1
    assert os.path.exists('results/perSamp/simSample2/infref_strain_dup.gff') == 1
    assert os.path.exists('results/perSamp/simSample2/infref_strain.fasta') == 1
    assert os.path.exists('results/perSamp/simSample2/infref_bowtie2_index') == 1
    assert os.path.exists('results/perSamp/simSample2/infref_strain.gff') == 1

    #Blast files
    assert os.path.exists('results/blast/blast.out') == 1
    assert os.path.exists('results/perSamp_blast/simSample1_blast.out') == 1
    assert os.path.exists('results/perSamp_blast/simSample2_blast.out') == 1

    #Trinity folders
    assert os.path.exists('results/perSamp_trinity/simSample1/') == 1
    assert os.path.exists('results/perSamp_trinity/simSample2/') == 1
    assert os.path.exists('results/trinity') == 1

    #Stat files
    assert os.path.exists('results/stats/samples.mapping.flagstat') == 1	
    assert os.path.exists('results/stats/simSample1.sorted.bam.flagstat') == 1	
    assert os.path.exists('results/stats/simSample2.sorted.bam.flagstat') == 1	

    #CDS coverage, gene coverage, depth and count files
    assert os.path.exists('results/coverage/infref/infref_genome_CDS_coverage.gct') == 1
    assert os.path.exists('results/coverage/infref/infref_genome_count.tsv') == 1
    assert os.path.exists('results/coverage/infref/infref_genome_depth.tsv') == 1
    assert os.path.exists('results/coverage/infref/infref_genome_gene_coverage.gct') == 1
    assert os.path.exists('results/coverage/inpt/inpt_genome_CDS_coverage.gct') == 1
    assert os.path.exists('results/coverage/inpt/inpt_genome_count.tsv') == 1
    assert os.path.exists('results/coverage/inpt/inpt_genome_depth.tsv') == 1
    assert os.path.exists('results/coverage/inpt/inpt_genome_gene_coverage.gct') == 1
    assert os.path.exists('results/coverage/simSample1/infref_genome_gene_coverage.gct') == 1
    assert os.path.exists('results/coverage/simSample1/infref_genome_CDS_coverage.gct') == 1
    assert os.path.exists('results/coverage/infref/infref_genome_depth_mean.tsv') == 1
    assert os.path.exists('results/coverage/simSample2/infref_genome_CDS_coverage.gct') == 1
    assert os.path.exists('results/coverage/infref/infref_genome_gene_coverage.gct') == 1
    assert os.path.exists('results/coverage/infref/infref_genome_depth_mqc.png') == 1
    

    #Variant calling files
    assert os.path.exists('results/variant-calling/infref/infref_simSample1.vcf') == 1
    assert os.path.exists('results/variant-calling/infref/infref_simSample2.vcf') == 1
    assert os.path.exists('results/variant-calling/infref/infref_simSample1_cleaned_allelicprimitives.vcf') == 1
    assert os.path.exists('results/variant-calling/infref/infref_simSample2_cleaned_allelicprimitives.vcf') == 1
    assert os.path.exists('results/variant-calling/inpt/inpt_simSample1.vcf') == 1
    assert os.path.exists('results/variant-calling/inpt/inpt_simSample2.vcf') == 1
    assert os.path.exists('results/variant-calling/inpt/inpt_simSample2_cleaned_allelicprimitives.vcf') == 1
    assert os.path.exists('results/variant-calling/inpt/inpt_simSample2_cleaned_allelicprimitives.vcf') == 1
    assert os.path.exists('results/variant-calling/inpt/infref_simSample2_cleaned_allelicprimitives.vcf') == 1
    assert os.path.exists('results/variant-calling/inpt/iinfref_simSample2_cleaned_allelicprimitives.vcf') == 1
    assert os.path.exists('results/variant-calling/perSamp/simSample1/simSample1.vcf') == 1
    assert os.path.exists('results/variant-calling/perSamp/simSample1/simSample1_cleaned_allelicprimitives.vcf') == 1
    assert os.path.exists('results/variant-calling/perSamp/simSample2/simSample2.vcf') == 1
    assert os.path.exists('results/variant-calling/perSamp/simSample2/simSample2_cleaned_allelicprimitives.vcf') == 1

def test_snakemake_variant_calling_results():

    #Check whether the correct variations have been detected for  the smaples

    simSample1_var = common.test_cleanvcf('results/variant-calling/infref/infref_simSample1_cleaned_allelicprimitives.vcf')
    simSample2_var = common.test_cleanvcf('results/variant-calling/inpt/inpt_simSample2_cleaned_allelicprimitives.vcf')
    assert simSample1_var == ['100', '1000']
    assert simSample2_var == ['200', '2000']
