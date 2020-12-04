*HBVouroboros* with freebayes variant calling and a RNA-seq simulator
===

## Freebayes variant calling

This branch of HBVouroboros now contains an additional snakemake file: "rules/freebayes_vc.smk". This is included in the main snakemake file and preforms variant calling using the freebayes software. The reuslts are stored in "results/var_valling_results". If you are using a conda environment from the HBVouroboros  master branch, you need to modify it to incude freebayes.

## RNA-seq simulator

An RNA-seq paired-end read simulator is added to the repository in the folder "RNAsim2". The program outputs two fastq.tar.gz files corresponding to left and right hand reads, as well as a sample annotation file. These files aredeposited in "RNAsim2/output". In order to run HBVouroboros using the simualted samples, the config file has to be modified to point to the corresponding sample annotation file. 

### Running the RNA-seq simualtor

The program needs to be run from the "RNAsim2/bin" folder. an example call looks like: 
```bash
python RNAsim.py 'gnl|hbvnuc|GQ924620_FT00000_C-C' '45000' '95' '45' --mutate --mutpos "66 88"
```

The above command will generate 45000 paired-end read, with the fragment size of 95 and read size of 45, using the genome  "GQ924620_FT00000_C-C" as reference. The reference genome will contain random mutations at positions 66 and 88. For more details and information further functions check the code. 


