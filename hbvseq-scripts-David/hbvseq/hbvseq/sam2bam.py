"""Convert SAM files to sorted BAM files"""

from .settings import hbvseqSettings

import tempfile
import re

def sam2bam(samfile, bamfile,
            keepNonMappingReads=False,
            sortBy='coordinate'):
    """Convert SAM files to sorted BAM files

    Also generates indexing (BAI) and flagstats"""

    if sortBy in ('readname', 'read', 'reads', 'name'):
        sortByFlag = '-n'
    elif sortBy in ('coordinate', 'coordinates', 'pos', 'position', 'coordinate'):
        sortByFlag = ''
    else:
        raise ValueError('sortBy must be either "readname" or "pos"')

    tempbam = tempfile.mkstemp(suffix='.bam', dir=hbvseqSettings.TMPDIR)[1]
    bampat = re.compile('(\.bam)')
    bai = bampat.subn('.bam.bai', bamfile)[0]
    flagstat = bampat.subn('.bam.flagstat', bamfile)[0]

    viewFlag = '' if keepNonMappingReads else '-F 4'

    samtools_map = ('{samtools} view {viewflag} '
                    '{input} -b -o {output}').format(samtools=hbvseqSettings.SAMTOOLS,
                                                     viewflag=viewFlag,
                                                     input=samfile,
                                                     output=tempbam)
    samtools_sort = ('{samtools} sort -@8 -T{prefix} -m5G {sortBy} '
                     '-O bam -o {bamfile} {input}').format(samtools=hbvseqSettings.SAMTOOLS,
                                                           prefix=tempbam,
                                                           sortBy=sortByFlag,
                                                           bamfile=bamfile,
                                                           input=tempbam)
    samtools_index = '{samtools} index {bamfile}'.format(samtools=hbvseqSettings.SAMTOOLS,
                                                         bamfile=bamfile)
    samtools_flagstat = ('{samtools} flagstat {bamfile} '
                         '> {flagstat}').format(samtools=hbvseqSettings.SAMTOOLS,
                                                bamfile=bamfile,
                                                flagstat=flagstat)
    chmod = 'chmod a+r,g+rw '+' '.join([bamfile, flagstat, bai])
    clean_files = 'rm -f {0}'.format(tempbam)

    commands = (samtools_map,
                samtools_sort,
                samtools_index,
                samtools_flagstat,
                chmod,
                clean_files)
    return(commands)
