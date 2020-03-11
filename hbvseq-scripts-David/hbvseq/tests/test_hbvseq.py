import unittest
import os
import re
import tempfile

import hbvseq


SUPPORTED_SUFFIX = ('.fastq', '.fastq.gz', '.fq', '.fq.gz',
                    '_1.fastq', '_1.fastq.gz', '_1.fq', '_1.fq.gz',
                    '_1_sequence.txt', '_1_sequence.txt.gz',
                    '_R1.fastq.gz', '_R1_1.fastq.gz',
                    '_R1_001.fastq.gz',
                    '_2.fastq', '_2.fastq.gz', '_2.fq', '_2.fq.gz',
                    '_2_sequence.txt', '_2_sequence.txt.gz',
                    '_R2.fastq.gz', '_R2_2.fastq.gz',
                    '_R2_001.fastq.gz')
SUPPORTED_SUFFIX_PES = [None]*4+['1']*9+['2']*9

ALPHAFASTQS = ['dir/myfastq{0}'.format(x)
               for x in SUPPORTED_SUFFIX]
ALPHAFASTQS_BASES = ['myfastq']*len(SUPPORTED_SUFFIX)
ALPHAFASTQS_PES = SUPPORTED_SUFFIX_PES

NUMRANGE = range(0, 20)
NUMFASTQS = ['dir/{0}{1}'.format(x, y)
             for x in NUMRANGE
             for y in SUPPORTED_SUFFIX]

NUMFASTQS_BASES = ['{0}'.format(x)
                   for x in NUMRANGE
                   for i in range(len(SUPPORTED_SUFFIX))]
NUMFASTQS_PES = SUPPORTED_SUFFIX_PES*len(NUMRANGE)

TESTFASTQS = ALPHAFASTQS+NUMFASTQS
TESTFASTQS_BASES = ALPHAFASTQS_BASES+NUMFASTQS_BASES
TESTFASTQS_PES = ALPHAFASTQS_PES+NUMFASTQS_PES




class TestIsFastqFile(unittest.TestCase):
    def test(self):
        for testfile in TESTFASTQS:
            self.assertEqual(hbvseq.isFastqFile(testfile), True)



def sam2bamHardcode(sam2bamResult,
                    keepNonMappingReads,
                    sortByReadName):
    if not isinstance(keepNonMappingReads, bool):
        raise ValueError('keepNonMappingReads must be bool')
    if not isinstance(sortByReadName, bool):
        raise ValueError('keepNonMappingReads must be bool')

    expLen = 6 if keepNonMappingReads else 8
    if len(sam2bamResult[0].split()) != expLen:
        print(sam2bamResult[0].split())
        raise ValueError('Inconsistent length with keepNonMappingReads')

    tmpFileInd = 5 if keepNonMappingReads else 7
    inFileInd = 2 if keepNonMappingReads else 4
    sortBy = '-n' if sortByReadName else ''
    
    infile = sam2bamResult[0].split()[inFileInd]
    bamfile = sam2bamResult[2].split()[-1]
    flagstat = bamfile+'.flagstat'
    bai = bamfile+'.bai'
    tmpfile = sam2bamResult[0].split()[tmpFileInd]

    flag = '' if keepNonMappingReads else '-F 4'
    commView = ('{samtools} view {flag} {infile} -b '
                '-o {tmpfile}').format(samtools="samtools",
                                       flag=flag,
                                       infile=infile,
                                       tmpfile=tmpfile)
    commSort = ('{samtools} sort -@8 '
                '-T{tmpfile} -m5G {sortBy} '
                '-O bam -o {bamfile} '
                '{tmpfile}').format(samtools="samtools",
                                    bamfile=bamfile,
                                    sortBy=sortBy,
                                    tmpfile=tmpfile)
    commIndex = ('{samtools} index '
                 '{bamfile}').format(samtools="samtools",
                                     bamfile=bamfile)
    commFlagStat = ('{samtools} flagstat {bamfile} '
                    '> {flagstat}').format(samtools="samtools",
                                           bamfile=bamfile,
                                           flagstat=flagstat)
    chmod = ('chmod a+r,g+rw '
             '{bamfile} {flagstat} {bai}').format(bamfile=bamfile,
                                                  flagstat=flagstat,
                                                  bai=bai)
    rm = 'rm -f {0}'.format(tmpfile)
    comms = (commView, commSort, commIndex, commFlagStat, chmod, rm)
    return comms


class TestSam2bam(unittest.TestCase):
    def test(self):
        insam = 'anydir/aSampleFile.sam'
        outbam = 'anyotherdir/aSampleFile.bam'
        sam2bamComms = hbvseq.sam2bam(insam, outbam,
                                      keepNonMappingReads=False,
                                      sortBy='readname')
        hardCodedComms = sam2bamHardcode(sam2bamComms,
                                         keepNonMappingReads=False,
                                         sortByReadName=True)
        self.assertEqual(hardCodedComms, sam2bamComms)

        # test in case not removing unmapped reads
        sam2bamCommsAll = hbvseq.sam2bam(insam, outbam,
                                         keepNonMappingReads=True,
                                         sortBy='readname')
        hardCodedSam2bamCommsAll = sam2bamHardcode(sam2bamCommsAll,
                                                   keepNonMappingReads=True,
                                                   sortByReadName=True)
        self.assertEqual(hardCodedSam2bamCommsAll, sam2bamCommsAll)

        # test sorting by position
        sam2bamCommsPos = hbvseq.sam2bam(insam, outbam,
                                         keepNonMappingReads=False,
                                         sortBy='pos')
        hardCodedSam2bamCommsPos = sam2bamHardcode(sam2bamCommsPos,
                                                   keepNonMappingReads=False,
                                                   sortByReadName=False)
        self.assertEqual(hardCodedSam2bamCommsPos, sam2bamCommsPos)


##class TestGsnapMapper(unittest.TestCase):
##    def test(self):
##        hsGnMap = hbvseq.GsnapMapper(genome=hbvseq.HUMAN_DEFGENOME,
##                                         genomeLabel=hbvseq.HUMAN_DEFGENOME_LABEL,
##                                         outdir='humanGs')
##        self.assertEqual(hsGnMap.__str__(),
##                         ('Mapper\n'
##                          '-- aligner:'
##                          '/apps64/bidev/ngs/gmap/gmap-2014-10-22/bin/gsnap\n'
##                          '-- aligner options:'
##                          '--gunzip --batch=4 --nthreads=12 --format=sam\n'
##                          'options controlling the behavior of the mapper'
##                          '-- outdir:humanGs\n'
##                          '-- keepNonMappingReads:True\n'
##                          '-- keepSam:False\n'))
##        self.assertEqual(hsGnMap.alignerBasename(), 'gsnap')
##        self.assertEqual(hsGnMap._mapSEcommand(fq1='test1_1.fq.gz',
##                                               samfile='test1.sam'),
##                         ('/apps64/bidev/ngs/gmap/gmap-2014-10-22/bin/gsnap '
##                          '-d humann_hg19_newIndex '
##                          '--gunzip --batch=4 --nthreads=12 --format=sam '
##                          'test1_1.fq.gz > humanGs/test1.sam'))
##        self.assertEqual(hsGnMap._mapPEcommand(fq1='test1_1.fq.gz',
##                                               fq2='test1_2.fq.gz',
##                                               samfile='test1.sam'),
##                         ('/apps64/bidev/ngs/gmap/gmap-2014-10-22/bin/gsnap '
##                          '-d humann_hg19_newIndex '
##                          '--gunzip --batch=4 --nthreads=12 --format=sam '
##                          'test1_1.fq.gz test1_2.fq.gz > humanGs/test1.sam'))
##        sampleSE = hbvseq.SeqSample('mysample', 'test0_1.fastq.gz', None)
##        samplePE = hbvseq.SeqSample('mysample',
##                                        'test0_1.fastq.gz', 'test0_2.fastq.gz')
##        self.assertEqual(hsGnMap.samFileName(sampleSE),
##                         'humanGs/mysample-hg19-gsnap.sam')
##        self.assertEqual(hsGnMap.samFileName(samplePE),
##                         'humanGs/mysample-hg19-gsnap.sam')
##        self.assertEqual(hsGnMap.bamFileName(sampleSE),
##                         'humanGs/mysample-hg19-gsnap.bam')
##        self.assertEqual(hsGnMap.bamFileName(samplePE),
##                         'humanGs/mysample-hg19-gsnap.bam')
##        hsGnMap.outdir = 'humanbam'
##        self.assertEqual(hsGnMap.map2samCommand(sampleSE),
##                         ('/apps64/bidev/ngs/gmap/gmap-2014-10-22/bin/gsnap '
##                          '-d humann_hg19_newIndex '
##                          '--gunzip --batch=4 --nthreads=12 --format=sam '
##                          'test0_1.fastq.gz '
##                          '> humanbam/mysample-hg19-gsnap.sam'))
##        self.assertEqual(hsGnMap.map2samCommand(samplePE),
##                         ('/apps64/bidev/ngs/gmap/gmap-2014-10-22/bin/gsnap '
##                          '-d humann_hg19_newIndex '
##                          '--gunzip --batch=4 --nthreads=12 --format=sam '
##                          'test0_1.fastq.gz test0_2.fastq.gz '
##                          '> humanbam/mysample-hg19-gsnap.sam'))
##
##        """test map2bamCommand"""
##        map2bamSE = hsGnMap.map2bamCommand(sampleSE)
##        self.assertEqual(len(map2bamSE), 8)
##        self.assertEqual(map2bamSE[0],
##                         ('/apps64/bidev/ngs/gmap/gmap-2014-10-22/bin/gsnap '
##                          '-d humann_hg19_newIndex '
##                          '--gunzip --batch=4 --nthreads=12 --format=sam '
##                          'test0_1.fastq.gz > '
##                          'humanbam/mysample-hg19-gsnap.sam'))
##        map2bamSEhardcode = sam2bamHardcode(map2bamSE[1:7],
##                                            keepNonMappingReads=True,
##                                            sortByReadName=False)
##        self.assertEqual(map2bamSE[1:7], list(map2bamSEhardcode))
##        self.assertEqual(map2bamSE[7],
##                         'rm -f humanbam/mysample-hg19-gsnap.sam')
##
##        map2bamPE = hsGnMap.map2bamCommand(samplePE)
##        self.assertEqual(len(map2bamPE), 8)
##        self.assertEqual(map2bamPE[0],
##                         ('/apps64/bidev/ngs/gmap/gmap-2014-10-22/bin/gsnap '
##                          '-d humann_hg19_newIndex '
##                          '--gunzip --batch=4 --nthreads=12 --format=sam '
##                          'test0_1.fastq.gz test0_2.fastq.gz > '
##                          'humanbam/mysample-hg19-gsnap.sam'))
##        map2bamPEhardcode = sam2bamHardcode(map2bamPE[1:7],
##                                            keepNonMappingReads=True,
##                                            sortByReadName=False)
##        self.assertEqual(map2bamPE[1:7], list(map2bamPEhardcode))
##        self.assertEqual(map2bamPE[7],
##                         'rm -f humanbam/mysample-hg19-gsnap.sam')
##

class TestBowtie2Mapper(unittest.TestCase):
    def test(self):
        hbvMap = hbvseq.Bowtie2Mapper(genome='abc',
                                      genomeLabel='HBV_abc',
                                      keepNonMappingReads=False,
                                      outdir='outbam')
        self.assertEqual(hbvMap.__str__(),
                         ('Mapper\n'
                          '-- aligner:bowtie2\n'
                          '-- aligner options:'
                          '--no-mixed --no-discordant -p 8\n'
                          'options controlling the behavior of the mapper'
                          '-- outdir:outbam\n'
                          '-- keepNonMappingReads:False\n'
                          '-- keepSam:False\n'))
        self.assertEqual(hbvMap._mapPEcommand(fq1='test1_1.fq.gz',
                                              fq2='test1_2.fq.gz',
                                              samfile='test1.sam'),
                         ('bowtie2 '
                          '--no-mixed --no-discordant -p 8 '
                          '-x abc '
                          '-1 test1_1.fq.gz -2 test1_2.fq.gz '
                          '-S outbam/test1.sam'))
        self.assertEqual(hbvMap._mapSEcommand(fq1='test1_1.fq.gz',
                                              samfile='test1.sam'),
                         ('bowtie2 '
                          '--no-mixed --no-discordant -p 8 '
                          '-x abc '
                          '-U test1_1.fq.gz -S outbam/test1.sam'))
        sampleSE = hbvseq.SeqSample('mysample', 'test0_1.fastq.gz', None)
        samplePE = hbvseq.SeqSample('mysample',
                                        'test0_1.fastq.gz',
                                        'test0_2.fastq.gz')
        self.assertEqual(hbvMap.map2samCommand(sampleSE),
                         ('bowtie2 '
                          '--no-mixed --no-discordant -p 8 '
                          '-x abc '
                          '-U test0_1.fastq.gz '
                          '-S outbam/mysample-HBV_abc-bowtie2.sam'))
        self.assertEqual(hbvMap.map2samCommand(samplePE),
                         ('bowtie2 '
                          '--no-mixed --no-discordant -p 8 '
                          '-x abc '
                          '-1 test0_1.fastq.gz -2 test0_2.fastq.gz '
                          '-S outbam/mysample-HBV_abc-bowtie2.sam'))

        # test map2bamCommand
        map2bamSE = hbvMap.map2bamCommand(sampleSE)
        self.assertEqual(map2bamSE[1],
                         ('bowtie2 '
                          '--no-mixed --no-discordant -p 8 '
                          '-x abc -U test0_1.fastq.gz '
                          '-S outbam/mysample-HBV_abc-bowtie2.sam'))
##       map2bamSEhardcode = sam2bamHardcode(map2bamSE[1:7],
##                                           keepNonMappingReads=False,
##                                           sortByReadName=False)
##
##        self.assertEqual(map2bamSE[1:7], list(map2bamSEhardcode))
##        self.assertEqual(map2bamSE[7],
##                         'rm -f outbam/mysample-HBV_abc-bowtie2.sam')

        map2bamPE = hbvMap.map2bamCommand(samplePE)
        self.assertEqual(map2bamPE[1],
                          ('bowtie2 '
                          '--no-mixed --no-discordant -p 8 '
                          '-x abc -1 test0_1.fastq.gz -2 test0_2.fastq.gz '
                          '-S outbam/mysample-HBV_abc-bowtie2.sam'))
##        map2bamPEhardcode = sam2bamHardcode(map2bamPE[2:8],
##                                            keepNonMappingReads=False,
##                                            sortByReadName=False)
##        self.assertEqual(map2bamPE[1:7], list(map2bamPEhardcode))
##        self.assertEqual(map2bamPE[7],
##                         'rm -f outbam/mysample-HBV_abc-bowtie2.sam')


class TestSequence(unittest.TestCase):
    def test(self):
        seq = hbvseq.Sequence('myseq', 'atcgc')
        self.assertEqual(seq.toFASTA(), '>myseq\natcgc')
        seq.extend(3)
        self.assertEqual(seq.toFASTA(), '>myseq extended by 3 bases\natcgcatc')
        self.assertRaises(ValueError, seq.extend, 9)
        tmpfile = tempfile.mkstemp()[1]
        seq.writeFASTA(tmpfile)
        with open(tmpfile, 'r') as f:
            inseq = f.readlines()
        self.assertEqual('>myseq extended by 3 bases\natcgcatc',
                         ''.join(inseq))
        seq2 = hbvseq.readSingleSeqFASTA(tmpfile)
        seq2.name = seq.name
        seq2.seq = seq.seq
        os.remove(tmpfile)


class TestFetchSeqCommand(unittest.TestCase):
    def test(self):
        myid = 'EM:U95551'
        myfile = 'anydir/U95551.fasta'
        expRes = ('ml load EMBOSS; seqret EM:U95551 -stdout -auto -osformat fasta '
                  '> anydir/U95551.fasta')
        self.assertEqual(expRes,
                         hbvseq.fetchSeqCommand(myid, myfile))


class TestFastaIndexName(unittest.TestCase):
    def test(self):
        self.assertEqual(hbvseq.fastaIndexName('myFasta.fa'),
                         'myFasta')
        self.assertEqual(hbvseq.fastaIndexName('myFasta2.fasta'),
                         'myFasta2')
        self.assertEqual(hbvseq.fastaIndexName('myFasta2.FASTA'),
                         'myFasta2')
        self.assertEqual(hbvseq.fastaIndexName('anydir/myFasta3.test'),
                         'anydir/myFasta3.test')


class TestFastqBaseNameExclPE(unittest.TestCase):
    def test(self):
        for ind in range(0, len(TESTFASTQS)):
            expRes = (TESTFASTQS_BASES[ind],
                      TESTFASTQS_PES[ind])
            self.assertEqual(hbvseq.getFastqBaseAndPE(TESTFASTQS[ind]),
                             expRes)


class TestParseFastqFilename(unittest.TestCase):
    def test(self):
        for ind in range(0, len(TESTFASTQS)):
            expRes = (TESTFASTQS_BASES[ind],
                      TESTFASTQS_PES[ind])
            self.assertEqual(hbvseq.parseFastqFilename(TESTFASTQS[ind]),
                             expRes)
        # define own ways
        self.assertEqual(hbvseq.parseFastqFilename('abc_Forward.fastq.gz',
                                                       '([a-z]*)_', 'Reverse'),
                         ('abc', '1'))
        self.assertEqual(hbvseq.parseFastqFilename('abc_Reverse.fastq.gz',
                                                       '([a-z]*)_', 'Reverse'),
                         ('abc', '2'))


class TestBowtie2IndexCommand(unittest.TestCase):
    def test(self):
        myFasta = 'anydir/myGenome.fasta'
        myFastaIndexComm = hbvseq.bowtie2IndexCommand(myFasta)
        expRes = ('ml load Bowtie2; bowtie2-build '
                  'anydir/myGenome.fasta anydir/myGenome')
        self.assertEqual(myFastaIndexComm,
                         expRes)

        myFasta2 = 'anydir/myGenome.fa'
        myFastaIndexComm2 = hbvseq.bowtie2IndexCommand(myFasta2,
                                                           'myIndex')
        expRes2 = ('ml load Bowtie2; bowtie2-build '
                   'anydir/myGenome.fa myIndex')
        self.assertEqual(myFastaIndexComm2,
                         expRes2)


class TestExtendSeq(unittest.TestCase):
    def test(self):
        seq = hbvseq.Sequence('myseq', 'atcgc')
        tmpfile = tempfile.mkstemp()[1]
        tmpfile2 = tempfile.mkstemp()[1]
        try:
            seq.writeFASTA(tmpfile)
            hbvseq.extendSeq(tmpfile, 5, tmpfile2)
            with open(tmpfile2, 'r') as f:
                self.assertEqual(''.join(f.readlines()),
                                 '>myseq extended by 5 bases\natcgcatcgc')
        finally:
            os.remove(tmpfile)
            os.remove(tmpfile2)


class TestEmbossid2name(unittest.TestCase):
    def test(self):
        self.assertEqual(hbvseq.embossid2name('EM:U99991'),
                         'U99991')
        self.assertEqual(hbvseq.embossid2name('EM:U99991', prefix='HBV-'),
                         'HBV-U99991')
        self.assertEqual(hbvseq.embossid2name('sw:bace1_human'),
                         'bace1_human')


class TestSeqSample(unittest.TestCase):
    def test(self):
        peSample = hbvseq.SeqSample('testID', 'test_1.fastq.gz', 'test_2.fastq.gz')
        self.assertEqual(peSample.isPE(), True)
        self.assertEqual(peSample.printComm('--fq1', '--fq2'),
                         '--fq1 test_1.fastq.gz --fq2 test_2.fastq.gz')
        seSample = hbvseq.SeqSample('testID', 'test_1.fastq.gz', None)
        self.assertEqual(seSample.isPE(), False)
        self.assertEqual(seSample.printComm('--fq1', None),
                         '--fq1 test_1.fastq.gz')


class TestSeqDataset(unittest.TestCase):
    def test(self):
        nTestFile = 10
        pefiles = ['testfile-{0}_combined_R{1}.fastq.gz'.format(x, y)
                   for x in range(0, nTestFile)
                   for y in (1, 2)]
        expSampleIDs = ['testfile-{0}'.format(x)
                        for x in range(0, nTestFile)]
        mySeqDataset = hbvseq.SeqDataset(pefiles,
                                  '^([\w|-]*)_combined',
                                  '_R2', outdir='/data64/bi/tmp')
        self.assertEqual(mySeqDataset.count(), nTestFile)
        self.assertEqual(mySeqDataset.isPE(), True)
        self.assertEqual(set(mySeqDataset.sampleIDs()),
                         set(expSampleIDs))
        for sample in mySeqDataset.seqSamples:
            ind = re.match("testfile-([0-9]*)", sample.id).group(1)
            self.assertEqual(sample.isPE(), True)
            self.assertEqual(sample.fq1,
                             'testfile-{0}_combined_R1.fastq.gz'.format(ind))
            self.assertEqual(sample.fq2,
                             'testfile-{0}_combined_R2.fastq.gz'.format(ind))


if __name__ == '__main__':
    unittest.main()
