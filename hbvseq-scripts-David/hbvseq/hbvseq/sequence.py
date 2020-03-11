"""A data structure to hold and extend sequence"""

from .settings import hbvseqSettings

import os.path
import subprocess
import re


class Sequence(object):
    """A simple data structure to hold and extend sequence"""

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def extend(self, N=1):
        if N < 1:
            raise ValueError('Extending length N must be no less than 1')
        elif N > len(self.seq):
            raise ValueError(('Extending length N must be no larger than '
                              'sequence length {0}').format(len(self.seq)))
        else:
            nseq = self.seq+self.seq[0:N]
            self.seq = nseq
            self.name = "{0} extended by {1} bases".format(self.name, N)

    def toFASTA(self):
        return '>{0}\n{1}'.format(self.name, self.seq)

    def writeFASTA(self, file):
        with open(file, 'w') as f:
            f.write(self.toFASTA())

def readSingleSeqFASTA(file):
    """Read a FASTA file where is only only one sequence"""
    name = ''
    seq = ''
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                name = line.lstrip('>').strip()
            else:
                seq += line.strip()
    res = Sequence(name, seq)
    return res


def fetchSeqCommand(id, file):
    """fetch sequence with emboss into a FASTA file"""
    fetch = ('ml load EMBOSS; {seqret} {id} -stdout -auto -osformat fasta '
             '> {file}').format(seqret=hbvseqSettings.SEQRET,
                                id=id, 
                                file=file)
    return fetch

def extendSeq(inFastaFile, N, outFastaFile):
    """Read a sequence from FASTA file, extends it by N sequences, and output into another FASTA file"""
    mySeq = readSingleSeqFASTA(inFastaFile)
    mySeq.extend(N)
    mySeq.writeFASTA(outFastaFile)

def fastaIndexName(fastaFile):
    fastqre = re.compile('(?i)(\.fa(sta)?)')
    fastaInd = fastqre.search(fastaFile, re.IGNORECASE)
    if fastaInd is None:
        indexName = fastaFile
    else:
        indexName = fastaFile[0:fastaInd.start()]
    return indexName

def embossid2name(embossid, prefix='', suffix=''):
    return '{prefix}{id}{suffix}'.format(prefix=prefix,
                                         id=embossid.split(':')[-1],
                                         suffix=suffix)

def bowtie2IndexCommand(inFastaFile, indexName=None):
    if indexName is None:
        indexName = fastaIndexName(inFastaFile)
    return 'ml load Bowtie2; {bowtie2build} {fasta} {index}'.format(bowtie2build=hbvseqSettings.BOWTIE2BUILD,
                                                   fasta=inFastaFile,
                                                   index=indexName)


def buildExtSeq(embossid, N, outdir='.'):
    """Extract sequence from EMBOSS, extend the sequence, and build indexes using bowtie2build

    The extended index name is returned.
    """
    idfilebase = embossid2name(embossid)
    idfile = os.path.join(outdir,
                          idfilebase+'.fasta')
    idextfile = os.path.join(outdir,
                             idfilebase+'_extended.fasta')
    fetchComm = fetchSeqCommand(embossid, idfile)
    subprocess.call(fetchComm, shell=True)
    # check whether seqret run successfully
    if not os.path.isfile(idfile) or os.path.getsize(idfile) == 0:
        raise ValueError(('seqfetch did not return the sequence of {}.'
                         'Check the identifier!').format(embossid))
    extendSeq(idfile, N, idextfile)
    extIndName = fastaIndexName(idextfile)
    bowtie2BuildComm = bowtie2IndexCommand(idextfile, extIndName)
    subprocess.call(bowtie2BuildComm, shell=True)
    bowtie2suffix = ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2',
                     '.rev.1.bt2', '.rev.2.bt2')
    bowtie2expfiles = [extIndName+suff for suff in bowtie2suffix]
    for bowtie2expfile in bowtie2expfiles:
        if not os.path.isfile(bowtie2expfile):
            msg = 'bowtie2 failed to build index {}'.format(bowtie2expfile)
            raise ValueError(msg)

    return extIndName
