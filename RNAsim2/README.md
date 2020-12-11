RNAsim2: a RNA-seq simulator to validate the sensitivity and specificity of
HBVouroboros
===

The `RNAsim2` program was written to validate the sensivity and specificity of
HBVouroboros.

It outputs two fastq.tar.gz files corresponding to left and right-hand reads, as
well as a sample annotation file. These files are deposited in "RNAsim2/output".

In order to run HBVouroboros using the simualted samples, the config file has to
be modified to point to the corresponding sample annotation file.

We are now in the process of integrating RNAsim2 in a snakemake workflow.

## Running RNAsim2

The program needs to be run from the "RNAsim2/bin" folder. An example call looks
like:

```bash
python RNAsim.py 'gnl|hbvnuc|GQ924620_FT00000_C-C' '45000' '95' '45' --mutate --mutpos "66 88"
```

The above command will generate 45000 paired-end read, with the fragment size of
95 and read size of 45, using the genome  "GQ924620_FT00000_C-C" as reference.
The reference genome will contain random mutations at positions 66 and 88.
