import os
import pytest
import argparse
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from subprocess import Popen, PIPE, STDOUT

sys.path.insert(0, 'workflow/rules/')
import common

#include: "rules/common.py"

def test_snakemake_output_files():
    #Run the pipeline
    p= Popen('snakemake --cores 20 --latency-wait 900', shell=True, stdout= PIPE, stderr= STDOUT)
    pout= p.stdout.read()
    print(pout.decode('utf-8'))
    
    #Check exits code and other expected output            
    #assert 0 == p.returncode
    #Read bam files
    assert os.path.exists('results/bam/simSample1.sorted.bam') == 1
    assert os.path.exists('results/bam/simSample2.sorted.bam') == 1
