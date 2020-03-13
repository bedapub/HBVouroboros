from Bio import Entrez
from Bio import SeqIO
from BCBio import GFF

Entrez.email = "jitao_david.zhang@roche.com"
with Entrez.efetch(
    db="nucleotide", rettype="gb", retmode="text", id="KC774468"
) as handle:
    seq_record = SeqIO.read(handle, "gb")

SeqIO.write(seq_record, "output.gb", "gb")

gb_handle = open('output.gb', 'r')
gff_handle = open('output.gff', 'w')
GFF.write(SeqIO.parse(gb_handle, "genbank"), gff_handle)
gff_handle.close()

