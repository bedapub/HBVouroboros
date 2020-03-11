"""Objects to make alignments"""

from .settings import hbvseqSettings
from .sam2bam import sam2bam

import os.path

class Mapper(object):
    """Base object of Mapper, a general object to host relevant information of an aligner"""

    aligner = ''
    aligneropts = ''

    # flags that influence the behavior of the object
    outdir = '.'
    keepNonMappingReads = False
    keepSam = False

    def __init__(self):
        pass

    def __str__(self):
        return ('Mapper\n'
                '-- aligner:{aligner}\n'
                '-- aligner options:{opts}\n'
                'options controlling the behavior of the mapper'
                '-- outdir:{outdir}\n'
                '-- keepNonMappingReads:{nm}\n'
                '-- keepSam:{sam}\n'.format(aligner=self.aligner,
                                            opts=self.aligneropts,
                                            outdir=self.outdir,
                                            nm=self.keepNonMappingReads,
                                            sam=self.keepSam))

    def alignerBasename(self):
        comm = self.aligner.split(";")[-1].strip() ## allow "XXX; comm"
        basename = os.path.basename(comm)
        return(basename)

    def validSam(self, samfile):
        if not samfile.startswith(self.outdir):
            samfile = os.path.join(self.outdir, samfile)
        return samfile

    def _mapSEcommand(self, fq1, samfile):
        """Map single-end FASTQ file to sam file

        samfile must have relative path to outdir"""
        pass

    def _mapPEcommand(self, fq1, fq2, samfile):
        """Map paired-end FASTQ files to sam file.

        samfile must have relative path to outdir"""
        pass

    def samFileName(self, seqSample):
        return os.path.join(self.outdir,
                            '{0}.sam'.format(seqSample.id))

    def bamFileName(self, seqSample):
        return os.path.join(self.outdir,
                            '{0}.bam'.format(seqSample.id))

    def map2samCommand(self, seqSample):
        samfile = self.samFileName(seqSample)

        if seqSample.fq2 is None:
            fq2sam = self._mapSEcommand(seqSample.fq1, samfile)
        else:
            fq2sam = self._mapPEcommand(seqSample.fq1, seqSample.fq2, samfile)

        return(fq2sam)

    def sam2bamCommands(self, seqSample):
        if not isinstance(self.keepNonMappingReads, bool):
            raise ValueError('keepNonMappingReads must be '
                             'either True or False')
        samfile = self.samFileName(seqSample)
        bamfile = self.bamFileName(seqSample)
        return sam2bam(samfile, bamfile,
                       keepNonMappingReads=self.keepNonMappingReads)

    def map2bamCommand(self, seqSample):
        """map FASTQ files to the genome and write the results in BAM format"""
        mod = hbvseqSettings.LOADMODULES
        comm1 = self.map2samCommand(seqSample)
        comms = self.sam2bamCommands(seqSample)
        res = [mod] + [comm1]+list(comms)
        if not self.keepSam:
            samfile = self.samFileName(seqSample)
            res += ['rm -f {samfile}'.format(samfile=samfile)]
        return res


class GenomeMapper(Mapper):
    """GenomeMapper is a specialised mapper used to align sequences to genomes"""

    samFormat = '{id}-{genome}-{aligner}.sam'
    bamFormat = '{id}-{genome}-{aligner}.bam'

    def __init__(self, genome, genomeLabel=None,
                 keepNonMappingReads=True, outdir='.'):
        if genomeLabel is None:
            genomeLabel = genome
        self.genome = genome
        self.genomeLabel = genomeLabel
        self.keepNonMappingReads = keepNonMappingReads
        self.outdir = outdir

    def samFileName(self, seqSample):
        keys = {'id': seqSample.id,
                'genome': self.genomeLabel,
                'aligner': self.alignerBasename()}
        return os.path.join(self.outdir,
                            self.samFormat.format(**keys))

    def bamFileName(self, seqSample):
        keys = {'id': seqSample.id,
                'genome': self.genomeLabel,
                'aligner': self.alignerBasename()}
        return os.path.join(self.outdir,
                            self.bamFormat.format(**keys))


class GsnapMapper(GenomeMapper):
    """GsnapMapper is a GenomeMapper using GSNAP as aligner"""
    aligner = hbvseqSettings.GSNAP
    aligneropts = hbvseqSettings.GSNAP_OPTS

    def __init__(self, genome,
                 genomeLabel=None,
                 keepNonMappingReads=True,
                 outdir='.'):
        super(GsnapMapper, self).__init__(genome,
                                          genomeLabel,
                                          keepNonMappingReads,
                                          outdir)

    def _mapSEcommand(self, fq1, samfile):
        fmt = '{prog} -d {genome} {args} {fq1} > {sam}'
        comm = fmt.format(prog=self.aligner,
                          args=self.aligneropts,
                          genome=self.genome,
                          fq1=fq1,
                          sam=self.validSam(samfile))
        return comm

    def _mapPEcommand(self, fq1, fq2, samfile):
        fmt = '{prog} -d {genome} {args} {fq1} {fq2} > {sam}'
        comm = fmt.format(prog=self.aligner,
                          args=self.aligneropts,
                          genome=self.genome,
                          fq1=fq1,
                          fq2=fq2,
                          sam=self.validSam(samfile))
        return comm


class Bowtie2Mapper(GenomeMapper):
    aligner = hbvseqSettings.BOWTIE2
    aligneropts = hbvseqSettings.BOWTIE2_OPTS

    def __init__(self, genome,
                 genomeLabel=None,
                 keepNonMappingReads=True,
                 outdir='.'):
        super(Bowtie2Mapper, self).__init__(genome,
                                            genomeLabel,
                                            keepNonMappingReads,
                                            outdir)

    def _mapSEcommand(self, fq1, samfile):
        """map SE file to a genome and write in SAM."""
        fmt = '{prog} {args} -x {genome} -U {fq1} -S {sam}'
        comm = fmt.format(prog=self.aligner,
                          args=self.aligneropts,
                          genome=self.genome,
                          fq1=fq1,
                          sam=self.validSam(samfile))
        return comm

    def _mapPEcommand(self, fq1, fq2, samfile):
        """map SE file to a genome and write in SAM."""
        fmt = '{prog} {args} -x {genome} -1 {fq1} -2 {fq2} -S {sam}'
        comm = fmt.format(prog=self.aligner,
                          args=self.aligneropts,
                          genome=self.genome,
                          fq1=fq1,
                          fq2=fq2,
                          sam=self.validSam(samfile))
        return comm
