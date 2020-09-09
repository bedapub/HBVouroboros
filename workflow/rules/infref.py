#!/usr/bin/env python3

import re
import os.path
from Bio import SeqIO
from Bio import Entrez
from BCBio import GFF


def get_infref_acc(blast_tab_file):
    """Get accession of the inferred reference strain

    Args:
        blast_tab_file (str): Tabular output file of BLAST (outfmt=6)

    Returns:
        str: accession number
    """

    with open(blast_tab_file, 'r') as f:
        return(f.read().split('\t')[1])

def get_infref_gb_acc(blast_tab_file):
    """Get GenBank accession of the inferred reference strain

    Args:
        blast_tab_file (str): Tabular output file of BLAST (outfmt=6)

    Returns:
        str: genbank accession number
    """

    acc = get_infref_acc(blast_tab_file)
    res = acc.split("|")[2].split("_")[0]
    return(res)

def download_gb(acc, outfile):
    """Download GenBank file with the given accession number

    Args:
        acc (str): GenBank acession number, example:'KC774468'
        outfile (str): output GenBank file

    Returns:
        int: number of records written as an integer.
    """

    res = -1
    Entrez.email = "jitao_david.zhang@roche.com"
    with Entrez.efetch(
      db="nuccore", rettype="gb", retmode="text", id=acc
    ) as handle:
      seq_record = SeqIO.read(handle, "gb")
    
    res = SeqIO.write(seq_record, outfile, 'gb')
    return(res)



def write_seq_by_acc(infile, acc, outfile):
    """Fetch sequence by accession number and write it to FASTA file

    Args:
        infile (str): FASTA file of genomes
        acc (str): Accession number of the inferred reference strain
        outfile (str): output FASTA file
    Returns:
        bool: if found or not
    """

    found = False
    seqs = SeqIO.parse(infile, 'fasta')
    for seq in seqs:
        if seq.id == acc:
            found = True
            SeqIO.write(seq, outfile, 'fasta')

    return(found)

def gb2gff(infile, outfile):
    """Translate GenBank file to GFF3 file. TODO: the procedure now does not 
    handle join correctly

    Args:
        infile (str): input GenBank file
        outfile (str): output GFF3 file
    Returns:
        Number of records written
    """

    gb_handle = open(infile, 'r')
    gff_handle = open(outfile, 'w')
    res = GFF.write(SeqIO.parse(gb_handle, "gb"), gff_handle)
    gff_handle.close()
    return(res)

def sort_FASTA_by_length(infile, outfile):
    """Sort FASTA sequence by length descendingly

    Args:
        infile (str): FASTA file
        outfile (str): Output FASTA file
    return:
        None
    """

    records = list(SeqIO.parse(infile, "fasta"))
    records.sort(key=lambda r: -len(r))
    SeqIO.write(records, outfile, "fasta")

def first_accession(fastafile):
    """Get accession of the first record in FASTA file

    Args:
        fastafile (str): FASTA file
    Returns:
        str : accession number of the first record
    """

    res = '';
    seqs=SeqIO.parse(fastafile, 'fasta')
    for seq in seqs:
        res = seq.id
        break;

    return(res)

def dup_gff(dup_fasta, ingff, outgff):
    """Make GFF files for duplicated genome

    Args:
        dup_fasta (str): duplicated FASTA
        ingff (str): input GFF file
        outgff (str): output GFF file
    Returns:
        None
    """

    acc = first_accession(dup_fasta)
    out_handle = open(outgff, 'w')
    with open(ingff, 'r') as f:
        for line in f:
            if line[0]=='#':
                out_handle.write(line)
            else:
                fl = line.split('\t')
                fl[0] = acc
                out_handle.write('\t'.join(fl))
            out_handle.write('\n')

    out_handle.close()
    return(None)
