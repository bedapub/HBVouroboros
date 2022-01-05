.. _Installation:

Installation
------------

Download the source code
########################

.. code-block:: console

   git clone https://github.com/bedapub/HBVouroboros.git


Setup conda environment
#######################

.. code-block:: console

   # setup conda environment
   cd envs; conda env create; cd -
   # in case it has been installed, use the command below to update
   # conda env update
   conda activate HBVouroboros


Run an example
##############

Test your installation by invoking pytest, in a terminal in the root directory of HBVouroboros:

.. code-block:: console

   $ pytest

This will run the pipeline with provided simulated RNA reads and test that the required output files are produced and mutations are correctly called.

