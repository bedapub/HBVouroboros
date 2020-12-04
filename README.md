*HBVouroboros* with freebayes variant calling and a RNA-seq simulator
===

## Freebayes variant calling

This branch of HBVouroboros now contains a further an additional snakemake file: "rules/freebayes_vc.smk". This is include in the main snakemake file and prefrms variant calling using the freebayes software. The reuslts are stored in "results/var_valling_results".

## RNAseq simulator

An RNA-seq paired-end read simulator is added to the repository inthe folder "RNAsim2". The program outputs two fastq.tar.gz files corresponding to left and right hand reads, as well as a sample annotation files. These files are to be found in "RNAsim2/output". In order to run HBVouroboros using the simualted samples, the config files has to be modified to point to the corresponding sample annotation file. 

### Running the RNA-seq simualtor

The program needs to be run from the "RNAsim2/bin" folder using a cokmand similar to 
```bash
python RNAsim.py 'gnl|hbvnuc|GQ924620_FT00000_C-C' '45000' '95' '45' --mutate --mutpos "66 88"```
```

The above command will generate 45000 paired-end read, with the fragment size of 95 and read size of 45, using the genome  "GQ924620_FT00000_C-C". The reads will contain random mutations at positions 66 and 88. For more details and information further functions check the code. 


