*HBVouroboros* automates sequencing-based HBV genotyping and expression profiling
===

This version addresses the issues 10, 11, 14.

## Issue 10: vcf files should remove/overlay duplicated genome parts
The vcf outputs of variant calling vcf output files of the form "results/variant calling/{reference genome}/{reference gemone}\_{sample}.vcf" are filtered and written to "{reference gemone}\_{sample}\_cleaned.vcf". These are then used to generate the aggrigated vcf.

## Issue 11: Use prespecified reference strain for variant calling and reporting
All analysis preformed using the inferred reference strain is now prefomred in parallel using an input reference strain. the input reference is specified using the paramter "inputRef" in "config.yaml". The directories in "reults" are modified to contain output correponding to inferred reference (infref) and input reference (inpt). Currently the input referene is not optional, if not specified the piple exits with an error. 

## Issue 14: Make RNAsim2 part of the snakemake pipeline
RNAsim2 can now be run with the pipeline by speifying corresponding paramers in "config.yaml". These parameters are to be found under the section "RNA simulation paramterts". If the paramter "doSim" is set to True. The RNAsim is ran and the pipeline is ran using the output.

## Issue 12: Allow genotype inference on either the study level or the sample level
This issue is not resolved. The current problem can be reproduced by setting "doSim" to false and uncommenting the last line in the main SnakeFile ( this will cause rules generating the error to be run). The error is generated when trinity is run with a BAM file of a single sample.

## Issue 12: VarScan vs Freebayes
Currently variant calling on simulated dat using freebayes does not report the SNPs, while mutations are called correctly when VarScan is used. To reproduce thie observation the following steps have to be taken. Set "doSim" in "config.yaml" to True and run HBVouroboros. Navigate to "results/variant calling/infref". The vcf output 
of variant calling does not show any SNPs. Naviage to "HBVouroboro/". And run the following commands to generate variant calling output using varScan:

```bash
samtools mpileup -f "results/infref/infref_strain_dup.fasta"  "results/infref_bam/infref_simSample.sorted.bam > mzData.mpileup"
java -jar VarScan.v2.3.9.jar pileup2snp myData.mpileup > varScanResults.txt
```
The current paramters in "config.yaml" specify 5 mutations. Varscan reports 10 mutations (5 in the original and 5 in the duplicated region) while Freebayes does not report anything.
