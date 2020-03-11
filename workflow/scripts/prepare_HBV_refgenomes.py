#!/usr/bin/env python3

import requests
import re
from Bio import SeqIO

def download_HBV_refgenomes(filename="hbvdbr.fas"):
  """Download HBV reference genomes and save them into a file

  Args:
      filename (str): The output filename

  Returns:
      bool: whether the request was ok
  """

  url = 'https://hbvdb.lyon.inserm.fr/data/references/hbvdbr.fas'
  r = requests.get(url)
  rstr = r.content.decode("utf-8")
  with open(filename, 'w') as f:
      f.write(rstr)

  return(r.ok)

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

if __name__ == '__main__':
    download_HBV_refgenomes("test.fasta") 
    dup_and_conc("test.fasta", "test-dupconc.fasta")
