#!/usr/bin/env python3

import argparse
import sys
from os import makedirs
import os.path
import snakemake
import pkg_resources
import yaml

trim_config = pkg_resources.resource_filename('HBVouroboros',
    'trim_reads/config/config.yaml')
default_primer_file = pkg_resources.resource_filename('HBVouroboros',
    'trim_reads/config/primers.fasta')
trim_snakefile = pkg_resources.resource_filename('HBVouroboros',
    'trim_reads/Snakefile')
trim_clusterfile = pkg_resources.resource_filename('HBVouroboros',
    'config/cluster.json')

if not os.path.exists(default_primer_file):
    raise Exception('default_primer_file not found: ' + default_primer_file)
if not os.path.exists(trim_config):
    raise Exception('trim_config not found: ' + trim_snakefile)
if not os.path.exists(trim_snakefile):
    raise Exception('trim_snakefile not found: ' + trim_snakefile)
if not os.path.exists(trim_clusterfile):
    raise Exception('trim_clusterfile not found: ' + trim_clusterfile)

def main(args):
    """start snakemake with input parameters"""
    sample_file = os.path.realpath(args.sample_annotation_file)
    with open(trim_config) as file:
       configs = yaml.load(file, Loader=yaml.FullLoader)
    local = args.local

    if not local:
        cluster_comm = ('sbatch -t {cluster.time} -c {cluster.cpu} '
                       '-N {cluster.nodes} --mem={cluster.mem} '
                       '--ntasks-per-node={cluster.ntasks_per_node}')
        cluster_config = trim_clusterfile
    else:
        cluster_comm = None
        cluster_config = None

    illumina_clip_file = args.illumina_clip_file
    if illumina_clip_file is None:
       illumina_clip_file = default_primer_file

    illumina_clip_opts = args.illumina_clip_opts
    if illumina_clip_opts is None:
       illumina_clip_opts = configs['illumina_clip_opts']

    trimmomatic_steps = args.trimmomatic_steps
    if trimmomatic_steps is None:
       trimmomatic_steps = configs['trimmomatic_steps']

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
            'illumina_clip_file': illumina_clip_file,
            'illumina_clip_opts': illumina_clip_opts,
            'trimmomatic_steps': trimmomatic_steps
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
    parser.add_argument('-l', '--local',
        action='store_true',
        help='If given, the commands will be executed locally. '
             'Otherwise, they will be submitted to the cluster with sbatch')
    parser.add_argument('--illumina_clip_file',
        default=None,
        help='A FASTA file containing adapters, primer sequences, etc. '
             'See Trimmomatic manual (\'fastaWithAdaptersEtc\') for details. '
             'If not provided, a default file of HBV-plasmid PCR is used.')
    parser.add_argument('--illumina_clip_opts',
        default=None,
        help='Options in combination with the clip file. '
             'If not provided, default settings are used '
             '(:2:30:10:2:keepBothReads).')
    parser.add_argument('--trimmomatic_steps',
        default=None,
        help='Further steps. If not provided, default settings are used '
             '(LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35 CROP:140).')

    args = parser.parse_args()

    sys.exit(main(args))

