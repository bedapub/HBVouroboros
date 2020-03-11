#!/usr/bin/env python3

import re
from Bio import SeqIO

def dup_and_conc_fasta(infile="hbvdbr.fas", outfile="hbvdbr-dupconc.bas"):
  """Duplicate sequences in a FASTA file, concatenate the original 
     with the duplicates, and write the concatenated sequences

  Args:
      infile (str): The input filename, pointing to a FASTA file
      outfile (str): The output filename, overwritten if the file exists.

  Returns:
      int: number of sequences 
  """

  count = 0
  fout = open(outfile, 'w')

  for record in SeqIO.parse(infile, 'fasta'):
    fout.write('>' + record.description + '\n')
    sstr = str(record.seq)
    fout.write(sstr*2 + '\n')
    count += 1

  fout.close()
  return(count)
