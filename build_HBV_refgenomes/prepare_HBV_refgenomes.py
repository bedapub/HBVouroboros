#!/usr/bin/env python3

import re
import os.path
from Bio import SeqIO

def get_simplified_id(desc):
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

def simplify_id(record):
    """Simplify IDs of HBV reference genomes, 
    using the format of GENOTYPE|Accession.

    Args:
        record (Bio.SeqRecord): A SeqRecord object
    Returns:
        Bio.SeqRecord
    """
    
    new_id = get_simplified_id(record.description)
    new_desc = new_id + ' ' + record.description
    record.id = new_id
    record.description = new_desc
    return(record)

def simplify_id_FASTA(infile, outfile):
    """Simplify IDs of all sequences in a FASTA file

    Args:
        infile (str): input file name
        outfile (str): output file name
    Returns:
        int: number of sequences written
    """

    sequences = SeqIO.parse(infile, 'fasta')
    outseqs = (simplify_id(record) for record in sequences)
    res = SeqIO.write(outseqs, outfile, 'fasta')
    return(res)

def dup_and_conc(record):
    """Duplicate the sequence and concatenate the original and duplicated
    sequence, and append a text label to the id and the description

    Args:
        record (Bio.SeqRecord): A SeqRecord object
    Returns:
        Bio.SeqRecord
    """
    record.seq = record.seq*2
    old_desc = record.description
    new_desc = old_desc.replace("length=", "original length=")
    new_desc += ' Duplicated and concatenated (final length:{})'.format(len(record.seq))
    record.description = new_desc
    record.id += '|DupConc'
    return(record)

def dup_and_conc_FASTA(infile, outfile):
    """Duplicate sequences in a FASTA file, concatenate the original 
       with the duplicates, and write the concatenated sequences

    Args:
        infile (str): The input filename, pointing to a FASTA file
        outfile (str): The output filename, overwritten if the file exists.

    Returns:
        int: number of sequences 
    """

    sequences = SeqIO.parse(infile, 'fasta')
    outseqs = (dup_and_conc(record) for record in sequences)
    res = SeqIO.write(outseqs, outfile, 'fasta')
    return(res)

def refgenome_ids(infile, replace_pipe=True):
    """Return IDs of reference genomes

    Args:
        infile (str): The input FASTA file name
        replace_pipe (bool): Whether replace pipe (|) with underscores (_)
    Returns:
        list: ids
    """
    sequences = SeqIO.parse(infile, 'fasta')
    res = (seq.id for seq in sequences)
    if replace_pipe: 
      res = (s.replace('|', '_') for s in res)
    return(res)

def split_FASTA(infile, outdir=None, prefix=''):
    """Split sequences in a FASTA file into separate files
       output file name is given by the ids (with pipes replaced by underscore)

       Args:
           infile (str): The input FASTA file name
           outdir (str): The output directory. Default: input file folder
           prefix (str): Prefix to the output file name
       Returns:
           int: number of sequences
    """

    count=0
    if outdir is None:
        outdir = os.path.dirname(infile)
    sequences = SeqIO.parse(infile, 'fasta')
    for seq in sequences:
        outfile = os.path.join(outdir, prefix + seq.id.replace('|', '_') + '.fasta')
        SeqIO.write(seq, outfile, 'fasta')
        count += 1

