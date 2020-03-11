# fetch reference genomes and make a merged pseudo genome
import os
from hbvseq.hbvseq import fetchSeqCommand

comms = []
files = []
with open('hbv-refgenomes.txt', 'r') as f:
    for line in f:
        if not (line.startswith('#') or
                line.startswith('Genotype')):
            genotype, id = line.strip().split('\t')
            usa = 'EM:'+id
            file = '{g}_{id}.fasta'.format(g=genotype,
                                           id=id)
            comms.append(fetchSeqCommand(usa, file))
            files.append(file)


for comm in comms:
    os.system(comm)

with open('merged-HBV-refgenomes.fasta', 'w') as fout:
    for file in files:
        with open(file, 'r') as fin:
            for line in fin:
                fout.write(line)
