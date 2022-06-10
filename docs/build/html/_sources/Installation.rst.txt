.. _Installation:

Installation
------------

Download the source code
########################

.. code-block:: console

   git clone https://github.com/bedapub/HBVouroboros.git


Setup conda environment
#######################

.. code-block:: python

   ## setup conda environment
   cd envs; conda env create; cd -
   ## in case it has been installed, use the command below to update
   ## conda env update
   conda activate HBVouroboros


Run an example
##############

Test your installation by invoking pytest, in a terminal in the root directory of HBVouroboros:

.. code-block:: console

   $ snakemake -j 99 --use-envmodules ## use --use-conda if no R module is present

