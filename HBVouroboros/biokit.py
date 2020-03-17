#!/usr/bin/env python3

from os.path import join, exists, isdir, isfile, expanduser, realpath

def biokit_sample_annotation_filename(biokit_outdir):
    return(join(biokit_outdir, 'samples.txt'))

def biokit_unmapped_sample_annotation(biokit_outdir, outfile):
    """ Get sample annotation from a biokit output directory

    Args:
        biokit_outdir (str): An output directory of the biokit pipeline
        outfile (str): Output file name of sample annotation
    return:
        number of samples
    """

    biokit_outdir = realpath(expanduser(biokit_outdir))
    if not isdir(biokit_outdir):
       raise Exception(
            'Directory {} does not exist'.format(biokit_outdir))
    infile = biokit_sample_annotation_filename(biokit_outdir)
    if not exists(infile):
       raise Exception(
           'sample annotaiton file ({}) not found'.format(infile))

    fout = open(outfile, 'w')
    with open(infile, 'r') as fin:
      header = fin.readline()
      fout.write(header)
      for line in fin:
          lsplit = line.rstrip().split('\t')
          fastq_format = join(
              biokit_outdir,
              'unmapped',
              '.'.join([
                   lsplit[0],
                  'unmapped_mate{}', 
                  'gz']))
          f1 = fastq_format.format(1)
          f2 = fastq_format.format(2)
          if not isfile(f1):
              raise Exception('File {} does not exist'.format(f1))
          if not isfile(f2):
              raise Exception('File {} does not exist'.format(f2))
          lsplit[2] = f1
          lsplit[3] = f2
          newline = '\t'.join(lsplit)
          fout.write(newline)
          fout.write('\n')

    fout.close()
