"""Objects to hold dataset"""

from .io import createDir
from .io import detectWarningPath
from .io import isFastqFile
from .io import parseFastqFilename

from .settings import hbvseqSettings

from .sequence import embossid2name
from .sequence import buildExtSeq

from .mapper import Bowtie2Mapper

from .jobs import executeCommandsOnCluster

import gzip
import shutil
import os.path
import tempfile
from subprocess import call

class SeqSample:
    """A sequencing sample"""
    def __init__(self, id, fq1, fq2=None):
        self.id = id
        self.fq1 = fq1
        if fq2 is not None and fq2 != '':
            self.fq2 = fq2
        else:
            self.fq2 = None

    def isPE(self):
        """Is the sample sequenced with paired ends?"""
        return self.fq2 is not None

    def printComm(self, fq1prefix='', fq2prefix=''):
        """Print file names for command-line tools"""
        res = ' '.join([fq1prefix, self.fq1])
        if(self.isPE()):
            res = ' '.join([res, fq2prefix, self.fq2])
        return res


class SeqDataset(object):
    """Create a SeqDataset object by parsing input directory or file list"""

    def __init__(self, input, idpattern=None, fq2pattern=None, outdir='.'):
        self.fastqfiles = []
        self.seqSamples = []
        self.isPEflag = False
        self.extLen = 500

        self.outdir = outdir
        self.logdir = os.path.join(outdir, 'logs')
        self.logfile = os.path.join(self.logdir, 'output.log')
        self.errfile = os.path.join(self.logdir, 'stderr.log')

        self.hbvMapper = None

        detectWarningPath(outdir)
        try:
            createDir(outdir)
            createDir(self.logdir)
        except IOError:
            print('IOError in creating outdir')

        self.parseInput(input, idpattern, fq2pattern)

    def isPE(self):
        return self.isPEflag

    def readlen(self, fileInd=0):
        if fileInd < 0 or fileInd > len(self.fastqfiles):
            raise ValueError('fileInd must be between 0 and '
                             '{}'.format(len(self.fastqfiles)))

        with gzip.open(self.fastqfiles[fileInd], 'r') as f:
            f.readline()  # skip first line
            return len(f.readline())-1

    def parseInput(self, input, idpattern, fq2pattern):
        if isinstance(input, str):
            if os.path.isfile(input):
                self.parseSampleInfoFile(input)
                return
            elif os.path.isdir(input):
                self.parseDir(input)
            else:
                raise ValueError('"Input" as a string must be a file '
                                 'or a directory')
        elif isinstance(input, (list, tuple)):
            for file in input:
                self.fastqfiles.append(file)

        self.parseFileNames(idpattern, fq2pattern)

    def parseSampleInfoFile(self, file):
        """Parse a sample info file, whcih contains four fields per line: sample ID, group, fq1, fq2

        Trailing spaces and comments (starting with #) are tolerated
        """
        self.isPEflag = True
        with open(file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue;
                fields = line.split('#')[0].rstrip().split('\t')
                self.seqSamples.append(SeqSample(id=fields[0],
                                                 fq1=fields[2],
                                                 fq2=fields[3]))

    def parseDir(self, dir):
        """Parse a directory iteratively to look up for FASTQ files"""
        for root, dirs, files in os.walk(dir):
            for entry in files:
                if isFastqFile(entry):
                    self.fastqfiles.append(os.path.join(dir, entry))

    def parseFileNames(self, idpattern, fq2pattern):
        sampleDict = {}
        for fastq in self.fastqfiles:
            ident, mpe = parseFastqFilename(fastq, idpattern, fq2pattern)
            if mpe == '2':
                self.isPEflag = True

            if ident not in sampleDict:
                sampleDict[ident] = dict([(mpe, fastq)])
            else:
                if mpe in sampleDict[ident]:
                    msg = ('More than one files detected for sample '
                           '{0}\nFASTA file number {1}\n Existing file {2}'
                           '\nNew file {3}').format(ident,
                                                    mpe,
                                                    sampleDict[ident][mpe],
                                                    os.path.basename(fastq))
                    raise ValueError(msg)
                sampleDict[ident][mpe] = fastq

        expCount = 2 if self.isPEflag else 1
        for sampleKey, sampleVal in sampleDict.items():
            if len(sampleVal) != expCount:
                raise ValueError('Sample {0} has {1} files -'
                                 '{2} expected'.format(sampleKey,
                                                       len(sampleVal),
                                                       expCount))
            self.seqSamples.append(SeqSample(id=sampleKey,
                                             fq1=sampleVal['1'],
                                             fq2=sampleVal['2']))

    def count(self):
        return len(self.seqSamples)

    def sampleIDs(self):
        ids = []
        for sample in self.seqSamples:
            ids.append(sample.id)
        return(ids)

    def __str__(self):
        str = '# SeqDataset with {count} samples\n'.format(count=self.count())
        str += '\t'.join(['id', 'fq1', 'fq2'])
        str += '\n'
        for sample in self.seqSamples:
            str += '\t'.join([sample.id, sample.fq1, sample.fq2])
            str += '\n'
        return str


class HBVSeqDataset(SeqDataset):

    refHBV = None
    refHBVlabel = None

    def __init__(self, input, idpattern=None, fq2pattern=None, outdir='.'):
        super(HBVSeqDataset, self).__init__(input,
                                            idpattern,
                                            fq2pattern,
                                            outdir)

        self.hbvExtGenomeDir = os.path.join(outdir, 'HBV-genome')
        self.hbvBamDir = os.path.join(outdir,
                                      'HBV-bamfiles')


    def _hbvDir(self, embossid):
        return os.path.join(self.outdir,
                            embossid2name(embossid,
                                          prefix='HBV-', suffix='-genome'))

    def _hbvBamDir(self, embossid):
        return os.path.join(self.outdir,
                            embossid2name(embossid,
                                          prefix='HBV-', suffix='-bamfiles'))

    def setRefHBV(self, embossid, forceRebuildIndex=False):
        self.refHBV = embossid
        self.refHBVlabel = embossid2name(embossid)
        
        refHBVdir = self._hbvDir(embossid)
        hbvBAMdir = self._hbvBamDir(embossid)

        self.hbvBamDir = hbvBAMdir

        build = False
        if not os.path.isdir(refHBVdir):
            try:
                createDir(refHBVdir)
            except:
                print("Directory of reference HBV genome failed to create!")
            build = True
        elif forceRebuildIndex:
            build = True

        if build:
            index = buildExtSeq(embossid,
                                N=self.extLen,
                                outdir=refHBVdir)

        if not os.path.isdir(hbvBAMdir):
            try:
                # note that buildExtSeq returns the correct index name
                createDir(hbvBAMdir)
            except ValueError:
                print('Failed to create BAM directory')

        if build:
            self.hbvMapper = Bowtie2Mapper(genome=index,
                                           genomeLabel=self.refHBVlabel,
                                           keepNonMappingReads=False,
                                           outdir=self.hbvBamDir)
        return
    
    def denovoHBV(self):
        """Map reads to reference HBV genomes and apply de-novo assembly"""
        pass
    
    def mapHBVcommands(self):
        if self.hbvMapper is None:
            raise LookupError('hbvMapper is not set! use setRefHBV with forceRebuildIndex as True!')
        comms = []
        for seqSample in self.seqSamples:
            comm = self.hbvMapper.map2bamCommand(seqSample)
            comms.append(";".join(comm))
        return(comms)

    def mapHBV(self):
        comms = self.mapHBVcommands()
        subComms = executeCommandsOnCluster(comms, cpu=hbvseqSettings.nthreads)
        return subComms

    def countMappedPairs(self):
        counts={}
        for s in self.seqSamples:
            id = s.id
            flagfile = id + "-" + self.refHBVlabel + "-bowtie2.bam.flagstat"
            if not os.path.isfile(os.path.join(self.hbvBamDir, flagfile)):
                raise ValueError("flagstat file for %s not found" % id)
            with open(os.path.join(self.hbvBamDir, flagfile)) as fo:
                for line in fo:
                    if(line.find("properly paired")>0):
                        pairCount=int(int(line.split(" ")[0])/2)
                        counts[id] = pairCount
        return(counts)

    def writeMappedPairs(self):
        counts = self.countMappedPairs();
        ofile = os.path.join(self.outdir,
                             "hbvseq-output-pairedReadCounts.tsv")
        with open(ofile, 'w') as f:
            f.write('{0}\t{1}\n'.format("Sample", "ReadPairCount"))
            [f.write('{0}\t{1}\n'.format(key, value)) for key,value in counts.items()]
        return([counts, ofile])

    def mpileup(self):
        bamfiles=[]
        for s in self.seqSamples:
            id = s.id
            bamfile = id + "-" + self.refHBVlabel + "-bowtie2.bam"
            bamfileFull = os.path.join(self.hbvBamDir, bamfile)
            if not os.path.isfile(bamfileFull):
                raise ValueError("flagstat file for %s not found" % id)
            bamfiles.append(bamfileFull)
        bamTmp = tempfile.NamedTemporaryFile(dir=hbvseqSettings.TMPDIR)
        with open(bamTmp.name, 'w') as f:
            for line in bamfiles:
                f.write(line)
                f.write('\n')
        mpileupRawOut = os.path.join(self.outdir,
                                     "hbvseq-output-mpileup-raw.txt")
        commRaw = ('ml load SAMtools; %s/samtools-mpielup-depth.Rscript '
                   '-infile %s -outfile %s') % (hbvseqSettings.BINDIR,
                                                bamTmp.name, mpileupRawOut)
        call(commRaw, shell=True)
        mpileupFinal = os.path.join(self.outdir,
                                    "hbvseq-output-mpileup-final.txt")
        commFinal = ('%s/depth_extended_sequence.Rscript ' +
                     '-infile %s -N %d -header -outfile %s') % (hbvseqSettings.BINDIR,
                                                                mpileupRawOut,
                                                                self.extLen,
                                                                mpileupFinal)
        call(commFinal, shell=True)
        return(mpileupFinal)
