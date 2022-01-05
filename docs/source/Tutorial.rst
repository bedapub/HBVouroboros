.. _Tutorial:

Turotial
------------

To run HBVouroboros follow these steps:


1. Prepare the SampleAnnotation file
#################################

list your input in a file named ``sampleAnnotation``. The file should be tab-limited and contain at least four columns. These are sample id (#ID), group (GROUP) and paths to R1 and R2 fastq files (FASTQ1 and FASTQ2). The default ``sampleAnnotation`` file is located at ``RNAsim2/output/sampleAnnotation``.

2. Prepare the config.yaml file
###############################

HBVouroboros is a snakemake pipeline and as such employs a config file located at ``config/config.yaml``, see related `snakemake documentation <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html>`_. Below you see the simplest set-up of ``config.yaml``.



.. code-block:: python

   ## use either way to indicate input files: either sample_annotation, or biokit_outdir
   ## sample_annotation should be a tab-delimited file containing at least four columns (#ID, GROUP, FASTQ1, FASTQ2)
   sample_annotation: RNAsim2/output/sampleAnnotation
   biokit_outdir: None

   ## options for read trimming with Trimmomatic
   trim_reads: False
   illumina_clip_file: 'resources/trim_reads/primers.fasta'
   illumina_clip_opts: ':2:30:10:2:keepBothReads'
   # trimmomatic steps (excluding ILLUMINACLIP)
   trimmomatic_steps: 'LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35 CROP:140'
   
   ## Run modes
   doInputRef: False
   inputRef: 'gnl|hbvnuc|AB064313_FT00000_C-G'
   doPerSamp: False
   
   ## RNAsim
   doSim: False
   sample_annotation_sm: RNAsim2/output/sampleAnnotation
   genomeId: gnl|hbvnuc|GQ924620_FT00000_C-C
   sampNum: "165000"
   pairedEndDist: "250"
   readLen: "150"
   mt: -mt
   mp: -mp 
   mtPos: "100"

Here only the parameter ``sample_annotation`` is of importance, it should be set to the path of your ``sampleAnnotation``. 


3. Choose Run-modes
###################

The program can be run in three different modes: 1) **infref**: the reference genome is inferred from the aggregation of all samples, 2) **inpt**: samples are compared against an input genotype of HBV, and 3) **per_Samp**: the reference genome is inferred for each sample. By default the pipeline runs the **infref** mode. The user can choose to run the **inpt** and **per_samp** modes by setting the relevant parameters ``doInputRef``, ``doPerSamp`` and ``inputRef`` in the config file.
 


4. Run snakemake
################

Run the pipeline by typing the following in a terminal at the project's root directory.

.. code-block:: console

   $ snakemake --cores 1 --use-conda

This will execute the pipeline localy. Please look at `snakemake documentation <https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_ for further information on how to run a snakemake pipeline.

Running snakemake with ``--forceall`` option
############################################

MultiQC adjusts filenames when it finds previous MultiQC output. If files names are changed as a result, snakemake will fail. This can occue when it is run with --forceall option and previous MultiQC html report(s) are present. Removing previous output in the ``results`` folder will solve this issue.


