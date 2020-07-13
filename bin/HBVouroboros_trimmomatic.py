#!/usr/bin/env python3

import argparse
import sys
from os import makedirs
import os.path
import snakemake
import pkg_resources


trim_snakefile = pkg_resources.resource_filename('HBVouroboros', 
    'trim_reads/Snakefile')
trim_clusterfile = pkg_resources.resource_filename('HBVouroboros',  
    'align_reads/config/cluster.json')

if not os.path.exists(trim_snakefile):
    raise Exception('trim_snakefile not found')
if not os.path.exists(trim_clusterfile):
    raise Exception('trim_clusterfile not found')

def main(args):

    sample_file = os.path.realpath(args.sample_annotation_file)
   
    cluster_comm = ('sbatch -t {cluster.time} -c {cluster.cpu} ' 
                   '-N {cluster.nodes} --mem={cluster.mem} ' 
                   '--ntasks-per-node={cluster.ntasks_per_node}')
    cluster_config = trim_clusterfile

    outdir = args.outdir
    if os.path.exists(outdir):
        makedirs(outdir, mode=0o775, exist_ok=True)

    status = snakemake.snakemake(
        snakefile=trim_snakefile,
        cluster=cluster_comm,
        cluster_config=cluster_config,
        cores=128, nodes=128, local_cores=4,
        config={
            'sample_annotation': sample_file,
            },
        workdir=outdir,
        restart_times=3,
        printshellcmds=True)

    if status: # translate "success" into shell exit code of 0
        return 0
    return 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Run trimmomatic prior to running HBVouroboros to trim and crop Illumina (FASTQ) data as well as to remove adaptors')

    parser.add_argument('sample_annotation_file',
        help = 'A sample annotation file, a tab-delimited file with a header line. '
               'The first column contains sample IDs, the second contains sample groups, the third forward read files, and '
               'the fourth reverse read files. Currently, single-read files are not tested.')
    parser.add_argument('--outdir', '-o',
        default='./HBVouroboros_trimmed',
        help='output directory: if missing, a folder "HBVouroboros_trimmed" '
             ' will be created in the current directory')

    args = parser.parse_args()

    sys.exit(main(args))

