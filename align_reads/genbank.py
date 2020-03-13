from Bio import Entrez
from Bio import SeqIO

Entrez.email = "jitao_david.zhang@roche.com"
with Entrez.efetch(
    db="nucleotide", rettype="gb", retmode="text", id="KC774468"
) as handle:
    seq_record = SeqIO.read(handle, "gb")

SeqIO.write(seq_record, "output.gb", "gb")
SeqIO.write(seq_record, "output.gff3", "bed")

