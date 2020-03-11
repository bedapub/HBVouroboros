"""Settings for hbvseq"""

import os

class Settings:
    
    nthreads = 8
    
    moduleNames =  ('Bowtie2', 'SAMtools', 'EMBOSS')
    LOADMODULES = ';'.join(['ml load %s' % mod for mod in moduleNames])

    ## executable programs
    GSNAP = 'gsnap'
    BOWTIE2 = 'bowtie2'
    BOWTIE2BUILD = 'bowtie2-build'
    SAMTOOLS = 'samtools'
    SEQRET = 'seqret'

    ## aligner parameters
    GSNAP_OPTS = ('--gunzip --batch=4 '
                      '--nthreads={} --format=sam').format(nthreads)
    BOWTIE2_OPTS = '--no-mixed --no-discordant -p {}'.format(nthreads)

    ## other parameters
    TMPDIR = '/pstore/data/bi/apps/tmp/'

    WARNING_PATHS = ['/homebasel/', '/pstore/home']

    BINDIR = '/pstore/apps/bioinfo/hbvseq/bin/'
    LOGDIR = '/pstore/data/bi/log/'
    LOGFORMAT = '%(asctime)s %(message)s'
    LOGFILE = os.path.join(LOGDIR,
                           'hbvseq.log')
    LOGDATEFORMAT = '%Y-%m-%d %H:%M:%S'

    def __init__(self):
        pass

hbvseqSettings = Settings()
