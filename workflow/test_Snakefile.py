import os
import pytest
import argparse
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from subprocess import Popen, PIPE, STDOUT
import yaml
sys.path.insert(0, 'workflow/rules/')
import common

with open('config/config.yaml', 'r') as stream:
	config = yaml.full_load(stream)
doInputRef = config['doInputRef']
doPerSamp = config['doPerSamp']


#This test is set to pass with the default setting of the config/config.yaml file

def test_snakemake_output_files():

    #Run the pipeline
    #p= Popen("snakemake --cores 10 --use-conda --latency-wait 900", shell=True, stdout= PIPE, stderr= STDOUT)
    #pout= p.stdout.read()
    #print(pout.decode('utf-8'))
     
    #Optionally, get stdout and stderr
    #stdout, stderr= p.communicate()

    #Check exits code and other expected output            
    #assert 0 == p.returncode

    #Read bam files
    assert os.path.exists('results/bam/Sample1.sorted.bam') == 1
    assert os.path.exists('results/bam/Sample2.sorted.bam') == 1

    #BAM and related files
    assert os.path.exists('results/infref_bam/infref_Sample1.sorted.bam.stat') == 1
    assert os.path.exists('results/infref_bam/infref_Sample1.sorted.bam') == 1
    assert os.path.exists('results/infref_bam/infref_Sample1.corrected.sorted.bam') == 1
    assert os.path.exists('results/infref_bam/infref_Sample1.sorted.bam.bai') == 1
    assert os.path.exists('results/infref_bam/infref_Sample2.sorted.bam.stat') == 1
    assert os.path.exists('results/infref_bam/infref_Sample2.sorted.bam') == 1
    assert os.path.exists('results/infref_bam/infref_Sample2.corrected.sorted.bam') == 1
    assert os.path.exists('results/infref_bam/infref_Sample2.sorted.bam.bai') == 1
    if doInputRef == True:
        assert os.path.exists('results/inpt_bam/inpt_Sample1.sorted.bam.stat') == 1
        assert os.path.exists('results/inpt_bam/inpt_Sample1.sorted.bam') == 1
        assert os.path.exists('results/inpt_bam/inpt_Sample1.corrected.sorted.bam') == 1
        assert os.path.exists('results/inpt_bam/inpt_Sample1.sorted.bam.bai') == 1
        assert os.path.exists('results/inpt_bam/inpt_Sample2.sorted.bam.stat') == 1
        assert os.path.exists('results/inpt_bam/inpt_Sample2.sorted.bam') == 1
        assert os.path.exists('results/inpt_bam/inpt_Sample2.corrected.sorted.bam') == 1
        assert os.path.exists('results/inpt_bam/inpt_Sample2.sorted.bam.bai') == 1
    if doPerSamp == True:
        assert os.path.exists('results/perSamp/Sample1/bam/Sample1.sorted.bam.stat') == 1
        assert os.path.exists('results/perSamp/Sample1/bam/Sample1.sorted.bam') == 1
        assert os.path.exists('results/perSamp/Sample1/bam/Sample1.corrected.sorted.bam') == 1
        assert os.path.exists('results/perSamp/Sample1/bam/Sample1.sorted.bam.bai') == 1
        assert os.path.exists('results/perSamp/Sample2/bam/Sample2.sorted.bam.stat') == 1
        assert os.path.exists('results/perSamp/Sample2/bam/Sample2.sorted.bam') == 1
        assert os.path.exists('results/perSamp/Sample2/bam/Sample2.corrected.sorted.bam') == 1
        assert os.path.exists('results/perSamp/Sample2/bam/Sample2.sorted.bam.bai') == 1

    #Strain fasta, gff, gb and bowtie index files
    assert os.path.exists('results/infref/infref_strain.gb') == 1
    assert os.path.exists('results/infref/infref_strain_dup.fasta') == 1
    assert os.path.exists('results/infref/infref_strain_dup.gff') == 1
    assert os.path.exists('results/infref/infref_strain.fasta') == 1
    assert os.path.exists('results/infref/infref_bowtie2_index') == 1
    assert os.path.exists('results/infref/infref_strain.gff') == 1
    if doInputRef == True:    
        assert os.path.exists('results/inpt/inpt_strain.gb') == 1
        assert os.path.exists('results/inpt/inpt_strain_dup.fasta') == 1
        assert os.path.exists('results/inpt/inpt_strain_dup.gff') == 1
        assert os.path.exists('results/inpt/inpt_strain.fasta') == 1
        assert os.path.exists('results/inpt/inpt_bowtie2_index') == 1
        assert os.path.exists('results/inpt/inpt_strain.gff') == 1
    if doPerSamp == True:
        assert os.path.exists('results/perSamp/Sample1/infref_strain.gb') == 1
        assert os.path.exists('results/perSamp/Sample1/infref_strain_dup.fasta') == 1
        assert os.path.exists('results/perSamp/Sample1/infref_strain_dup.gff') == 1
        assert os.path.exists('results/perSamp/Sample1/infref_strain.fasta') == 1
        assert os.path.exists('results/perSamp/Sample1/infref_bowtie2_index') == 1
        assert os.path.exists('results/perSamp/Sample1/infref_strain.gff') == 1
        assert os.path.exists('results/perSamp/Sample2/infref_strain.gb') == 1
        assert os.path.exists('results/perSamp/Sample2/infref_strain_dup.fasta') == 1
        assert os.path.exists('results/perSamp/Sample2/infref_strain_dup.gff') == 1
        assert os.path.exists('results/perSamp/Sample2/infref_strain.fasta') == 1
        assert os.path.exists('results/perSamp/Sample2/infref_bowtie2_index') == 1
        assert os.path.exists('results/perSamp/Sample2/infref_strain.gff') == 1

    #Blast files
    assert os.path.exists('results/blast/blast.out') == 1
    if doInputRef == True: 
        assert os.path.exists('results/perSamp_blast/Sample1_blast.out') == 1
    if doPerSamp == True:
        assert os.path.exists('results/perSamp_blast/Sample2_blast.out') == 1

    #Trinity folders
    if doPerSamp == True:
        assert os.path.exists('results/perSamp_trinity/Sample1/') == 1    
        assert os.path.exists('results/perSamp_trinity/Sample2/') == 1
    
    assert os.path.exists('results/trinity') == 1

    #Stat files
    assert os.path.exists('results/stats/samples.mapping.flagstat') == 1
    assert os.path.exists('results/stats/Sample1.sorted.bam.flagstat') == 1	
    assert os.path.exists('results/stats/Sample2.sorted.bam.flagstat') == 1	

    #CDS coverage, gene coverage, depth and count files
    assert os.path.exists('results/coverage/infref/infref_genome_CDS_coverage.gct') == 1
    assert os.path.exists('results/coverage/infref/infref_genome_count.tsv') == 1
    assert os.path.exists('results/coverage/infref/infref_genome_depth.tsv') == 1
    assert os.path.exists('results/coverage/infref/infref_genome_gene_coverage.gct') == 1
    if doInputRef == True:        
        assert os.path.exists('results/coverage/inpt/inpt_genome_CDS_coverage.gct') == 1
        assert os.path.exists('results/coverage/inpt/inpt_genome_count.tsv') == 1
        assert os.path.exists('results/coverage/inpt/inpt_genome_depth.tsv') == 1
        assert os.path.exists('results/coverage/inpt/inpt_genome_gene_coverage.gct') == 1
    if doPerSamp == True:    
        assert os.path.exists('results/coverage/Sample1/genome_gene_coverage.gct') == 1
        assert os.path.exists('results/coverage/Sample1/genome_CDS_coverage.gct') == 1
        assert os.path.exists('results/coverage/Sample2/genome_gene_coverage.gct') == 1
        assert os.path.exists('results/coverage/Sample2/genome_CDS_coverage.gct') == 1       
        assert os.path.exists('results/coverage/Sample1/genome_count.tsv') == 1
        assert os.path.exists('results/coverage/Sample1/genome_depth.tsv') == 1
        assert os.path.exists('results/coverage/Sample2/genome_count.tsv') == 1
        assert os.path.exists('results/coverage/Sample2/genome_depth.tsv') == 1    
    assert os.path.exists('results/coverage/infref/infref_genome_depth_mqc.png') == 1
    assert os.path.exists('results/coverage/infref/infref_genome_depth_mean.tsv') == 1
    if doInputRef == True:        
        assert os.path.exists('results/coverage/inpt/inpt_genome_depth_mqc.png') == 1
        assert os.path.exists('results/coverage/inpt/inpt_genome_depth_mean.tsv') == 1

    #AA Variant calling files
    assert os.path.exists('results/variant-calling-AA/infref/infref_Sample1.sorted.sam') == 1    
    assert os.path.exists('results/variant-calling-AA/infref/infref_Sample2.sorted.sam') == 1
    if doInputRef == True:    
        assert os.path.exists('results/variant-calling-AA/inpt/inpt_Sample1.sorted.sam') == 1    
        assert os.path.exists('results/variant-calling-AA/inpt/inpt_Sample2.sorted.sam') == 1
    if doPerSamp == True:    
        assert os.path.exists('results/variant-calling-AA/perSamp/Sample1/Sample1.sorted.sam') == 1
        assert os.path.exists('results/variant-calling-AA/perSamp/Sample2/Sample2.sorted.sam') == 1

    #Variant calling files
    assert os.path.exists('results/variant-calling/infref/infref_Sample1_varscan.vcf') == 1
    assert os.path.exists('results/variant-calling/infref/infref_Sample2_varscan.vcf') == 1
    assert os.path.exists('results/variant-calling/infref/infref_Sample1_cleaned_allelicprimitives.vcf') == 1
    assert os.path.exists('results/variant-calling/infref/infref_Sample2_cleaned_allelicprimitives.vcf') == 1
    if doInputRef == True:        
        assert os.path.exists('results/variant-calling/inpt/inpt_Sample1_varscan.vcf') == 1
        assert os.path.exists('results/variant-calling/inpt/inpt_Sample2_varscan.vcf') == 1
        assert os.path.exists('results/variant-calling/inpt/inpt_Sample1_cleaned_allelicprimitives.vcf') == 1
        assert os.path.exists('results/variant-calling/inpt/inpt_Sample2_cleaned_allelicprimitives.vcf') == 1
    if doPerSamp == True:    
        assert os.path.exists('results/variant-calling/perSamp/Sample1/Sample1_varscan.vcf') == 1
        assert os.path.exists('results/variant-calling/perSamp/Sample1/Sample1_cleaned_allelicprimitives.vcf') == 1
        assert os.path.exists('results/variant-calling/perSamp/Sample2/Sample2_varscan.vcf') == 1
        assert os.path.exists('results/variant-calling/perSamp/Sample2/Sample2_cleaned_allelicprimitives.vcf') == 1



def test_snakemake_variant_calling_results():

    #Check whether the correct variations have been detected for  the smaples

    Sample1_var = common.test_cleanvcf('results/variant-calling/infref/infref_Sample1_cleaned_allelicprimitives.vcf','results/infref/infref_strain_dup.fasta')
    if doInputRef == True:     
        Sample2_var = common.test_cleanvcf('results/variant-calling/inpt/inpt_Sample2_cleaned_allelicprimitives.vcf', 'results/inpt/inpt_strain_dup.fasta')
    assert Sample1_var == ['925', '934', '1371', '1896']
    #assert simSample2_var == ['582', '925', '928', '934', ]
