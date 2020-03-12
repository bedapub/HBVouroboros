#!/usr/bin/env python3

import re
from Bio import SeqIO

def simplified_id(desc):
   """ Make a new id for reference genome that contains genotype and accession

   Args:
      desc (str): The input description
   
   Returns:
      str: a new id
   """
   
   strain = desc.split(" ")[0].split("|")[2]
   acc = strain.split("_")[0]
   gt = strain.split("-")[1]
   id = "{}|{}".format(gt, acc)

   return(id)

def dup_and_conc(infile="hbvdbr.fas", outfile="hbvdbr-dupconc.bas"):
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
    newid = simplified_id(record.description)
    fout.write('>' + newid + ' ' + record.description + ' Duplicated_Concatenated\n')
    sstr = str(record.seq)
    fout.write(sstr*2 + '\n')
    count += 1

  fout.close()
  return(count)
