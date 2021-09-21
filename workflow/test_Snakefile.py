import os
import pytest
import argparse
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess as sp

sys.path.insert(0, 'workflow/rules/')
import common

#include: "rules/common.py"

def test_snakemake_output_files():
    #Run the pipeline
    p= sp.Popen('snakemake --cores 4 --forceall', shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
     
    #Optionally, get stdout and stderr
    stdout, stderr= p.communicate()

    #Check exits code and other expected output            
    #assert 0 == p.returncode
