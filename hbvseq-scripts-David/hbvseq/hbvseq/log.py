"""Logging facility"""
import logging
import os.path

from .settings import hbvseqSettings

logging.basicConfig(
    level=logging.DEBUG,
    format=hbvseqSettings.LOGFORMAT,
    filename=hbvseqSettings.LOGFILE
)


def startLog(filename, loglevel='INFO'):
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(filename=os.path.join(LOGDIR, filename),
                        level=numeric_level,
                        format=hbvseqSettings.LOGFORMAT, 
                        datefmt=hbvseqSettings.LOGDATEFORMAT)
    logging.info('Program started by {0}'.format(os.environ['LOGNAME']))


def endLog():
    logging.info('Program exits')
