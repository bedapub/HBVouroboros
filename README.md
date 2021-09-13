*HBVouroboros* automates sequencing-based HBV genotyping and expression profiling
===

*HBVouroboros* uses RNA-sequencing reads to infer HBV genotype, quantify HBV
transcript expression, and perform variant calling of HBV genomes.

*HBVouroboros*, distributed under the GPL-3 license, is available at
https://github.com/bedapub/HBVouroboros.

## Installation and usage

### Download the source code

```bash
git clone https://github.com/bedapub/HBVouroboros.git
```

### Setup conda environment

```bash
## setup conda environment
cd envs; conda env create; cd -
## in case it has been installed, use the command below to update
## conda env update
conda activate HBVouroboros
```

### Run an example

An out-of-box example can be run by starting the `snakemake` pipeline.

```bash
snakemake -j 99 --use-envmodules ## use --use-conda if no R module is present
```

### Run the pipeline with your own data

Modify the `config/config.yaml` file to specify a sample annotation file.

### Run HBVouroboros using unmapped reads from a Biokit output directory

This feature has been disabled now. It may be activated in the future.

### Validating the sensitivity and specificity of HBVouroboros with RNAsim2

We created RNAsim2, a RNA-seq simulator to validate the sensitivity and
specificity of HBVouroboros. See [RNAsim2/README.md](RNAsim2/README.md) for
details.

## New version 13_09_2-21


### Variant calling against a circular genome

While the duplication solution works well for reference genome inference, there are implications if it used in the sam way for variant calling purposes. A duplicated genome contrains two mapping possibilities for reads not spanning over the either ends of the linear version of the genome. This influences the confidence in mapping of reads with a position-based bias. To improve this, we edit the bam files mapped against the duplicated genome by shifting back the reads whose leftmost position is in the duplicated region, to their corresponding location on the single genome. However after the shift, the reads that span the end of the genome will now map to locations beyond the original length. Therefore the reads are mapped to a genome that is longer that the orignal by up to the maximum length of the reads. This means reads within the range of the maximum read length from the origin will map to two locations.Tthe outcome is a relatively lower precision in variant calling in the aformentioned region of the genome which would be only noticable with significantly low read coverage.

### Variant calling using an input reference

The pipeline can be run with a user-provided reference genome. To do so, in the config file set the parameter 'doInputRef' to 'True' and set 'inputRef' to a genome present in the list of reference genomes provided in the folder 'resources'.

### inference of reference gneome on a per sample basis

The pipeline can be used to infere reference gneomes at the level of each sample. To do so, in the config file set the parameter 'doPerSamp' to true.


### Simulated samples for testing the new version

The folder 'simSamples' is added th the 'RNAsim2' folder. It contains two simualted samples for testing the pipline using the default parameters in the config file. The samples are as generated using these commands:

* python RNAsim.py 'gnl|hbvnuc|GQ924620_FT00000_C-C' '90000' '95' '45' --mutate --mutpos "100 1000"
* python RNAsim.py 'gnl|hbvnuc|AB064313_FT00000_C-G' '90000' '95' '45' --mutate --mutpos "200 2000"

For informaiton on RNAsim parameters please see the README.md file in the 'RNAsim2' folder.



## Known issues and solutions

### What to do if conda environment initialization takes too long?

Above we use the default conda solver. If you suffer from slow speed of conda,
consider using [mamba](https://github.com/mamba-org/mamba), which is a drop-in
replacement of conda.

If you met more issues, please raise them using the Issues function of GitHub.

