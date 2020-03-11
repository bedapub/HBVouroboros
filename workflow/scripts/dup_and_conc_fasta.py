#!/usr/bin/env python3

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
  with open(outfile, 'w') as fout:
    with open(infile, 'r') as fin:
      while True:
        line = fin.readline().strip()
        
        if len(line)==0 or not line:
            break;

        if line[0]=='>':
          count += 1
          fout.write(line)
          fout.write('\n')
        else:
          conc = line + line + '\n'
          fout.write(conc)

  return(count)
