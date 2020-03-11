"""I/O functions to parse input files"""

from .settings import hbvseqSettings

import warnings
import os
import errno
import re

FASTQPAT = re.compile('(?:_|R|(?:_R))?((?<=[_|R])[1|2])?'
                      '((?<=R[1|2])_0{0,2}[1|2])?'
                      '\.f(ast)?q(\.gz)?$')
SEQTXTPAT = re.compile('_?([1|2])?_sequence.txt(\.gz)?$')


def createDir(dir):
    """Create directory if not exists, and raise error if the creation fails"""
    try:
        os.makedirs(dir, mode=0o775)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise



def detectWarningPath(dir):
    """Detect paths where output should not be saved"""
    absdir = os.path.abspath(dir)
    for warnPath in hbvseqSettings.WARNING_PATHS:
        if re.search(warnPath, absdir) is not None:
            warnMessage = ('Directory {0} found matching the warning path {1}.'
                           'Are you sure about this?').format(dir, warnPath)
            warnings.warn(warnMessage, UserWarning)


def isFastqFile(filename):
    """Test whether a file name looks like FASTQ file"""
    return (not filename.startswith('.') and
            (FASTQPAT.search(filename) is not None or
             SEQTXTPAT.search(filename) is not None))


def getFastqBaseAndPE(filename):
    """
    Get basename and PE information from FASTQ file names
    """
    filebase = os.path.basename(filename)

    txtSearch = SEQTXTPAT.search(filebase)
    fastqSearch = FASTQPAT.search(filebase)

    if txtSearch is not None:
        base = filebase[0:txtSearch.start()]
        pe = txtSearch.group(1)
    elif fastqSearch is not None:
        base = filebase[0:fastqSearch.start()]
        pe = fastqSearch.group(1)
    else:
        raise ValueError('Bug in getFastqBaseAndPE detected.'
                         'File name:{0}'.format(filename))
    return base, pe


def parseFastqFilename(fastqFileName, idPatStr=None, fq2PatStr=None):
    """Get sample identifier and PE information from a FASTQ file name"""
    fastqbase = os.path.basename(fastqFileName)
    guess = getFastqBaseAndPE(fastqbase)
    if idPatStr is None and fq2PatStr is None:
        return guess
    guess = list(guess)
    if idPatStr is not None:
        idmatch = re.search(idPatStr, fastqbase)
        if idmatch is None:
            warnings.warn('id match failed: id={0}; file={1}.'
                             'Use best guess instead'.
                             format(idPatStr, fastqbase))
        else:
            guess[0] = ''.join(idmatch.groups())
    if fq2PatStr is not None:
        fq2match = re.search(fq2PatStr, fastqbase)
        if fq2match is not None:
            guess[1] = '2'
        else:
            guess[1] = '1'
    return tuple(guess)
