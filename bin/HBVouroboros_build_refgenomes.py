#!/usr/bin/env python3

import sys
from os import makedirs
import os.path
import snakemake
import pkg_resources
import argparse

refgenome_snakefile = pkg_resources.resource_filename('HBVouroboros', 
    'build_refgenomes/Snakefile')

def main(args):
    outdir = args.outdir
    if outdir is None:
         makedirs(outdir, mode=0x775, exist_ok=True)

    status = snakemake.snakemake(
        snakefile = refgenome_snakefile,
        cores=2, nodes=2, local_cores=2,
        workdir=outdir,
        printshellcmds=True)

    if status:
        return 0
    return 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Download HBV genomes and build Bowtie2 indices')

    parser.add_argument('outdir',
        help = 'Output directory of reference genomes')

    args = parser.parse_args()
    sys.exit(main(args))
