#!/usr/bin/env python3

import argparse
import sys
from os import makedirs
import os.path
import snakemake
import pkg_resources
import yaml

pisces_config = pkg_resources.resource_filename('HBVouroboros',
    'pisces/config/config.yaml')
pisces_snakefile = pkg_resources.resource_filename('HBVouroboros',
    'pisces/Snakefile')
clusterfile = pkg_resources.resource_filename('HBVouroboros',
    'config/cluster.json')

if not os.path.exists(pisces_config):
    raise Exception('pisces_config not found: ' + pisces_config)
if not os.path.exists(pisces_snakefile):
    raise Exception('pisces_snakefile not found: ' + pisces_snakefile)
if not os.path.exists(clusterfile):
    raise Exception('clusterfile not found: ' + clusterfile)

def main(args):
    """start snakemake with input parameters"""
    sample_file = os.path.realpath(args.sample_annotation_file)
    ref_genome_fa = os.path.realpath(args.ref_genome_fa)
    local = args.local

    if not local:
        cluster_comm = ('sbatch -t {cluster.time} -c {cluster.cpu} '
                       '-N {cluster.nodes} --mem={cluster.mem} '
                       '--ntasks-per-node={cluster.ntasks_per_node}')
        cluster_config = clusterfile
    else:
        cluster_comm = None
        cluster_config = None

    outdir = args.outdir
    if os.path.exists(outdir):
        makedirs(outdir, mode=0o775, exist_ok=True)

    status = snakemake.snakemake(
        snakefile=pisces_snakefile,
        cluster=cluster_comm,
        cluster_config=cluster_config,
        cores=128, nodes=128, local_cores=4,
        config={
            'sample_annotation': sample_file,
            'ref_genome_fa': ref_genome_fa,
            'create_genome_size_cmd': args.create_genome_size_cmd,
            'pisces_cmd': args.pisces_cmd,
            'pisces_opts': args.pisces_opts,
            'ref_genome_genus': args.genus,
            'ref_genome_species': args.species,
            'ref_genome_build': args.build
            },
        workdir=outdir,
        restart_times=3,
        printshellcmds=True)

    if status: # translate "success" into shell exit code of 0
        return 0
    return 1

if __name__ == '__main__':

    with open(pisces_config) as file:
       def_configs = yaml.load(file, Loader=yaml.FullLoader)

    parser = argparse.ArgumentParser(
        description = 'Run Pisces for variant calling.')

    parser.add_argument('sample_annotation_file',
        help = 'A sample annotation file, a tab-delimited file with a header line. '
               'The first column contains sample IDs, the second contains sample groups, the third forward read files, and '
               'the fourth reverse read files. Currently, single-read files are not tested.')
    parser.add_argument('ref_genome_fa',
        help = 'FASTA sequences of the reference genome')

    parser.add_argument('--outdir', '-o',
        default='./HBVouroboros_trimmed',
        help='output directory: if missing, a folder "HBVouroboros_trimmed" '
             ' will be created in the current directory')
    parser.add_argument('-l', '--local',
        action='store_true',
        help='If given, the commands will be executed locally. '
             'Otherwise, they will be submitted to the cluster with sbatch')
        
    parser.add_argument('--create-genome-size-cmd',
        default=def_configs['create_genome_size_cmd'],
        help='The command to create the GenomeSize file for Pisces.')
    parser.add_argument('--pisces-cmd',
        default=def_configs['pisces_cmd'],
        help='The command to run the Pisces program.')
    parser.add_argument('--pisces-opts',
        default=def_configs['pisces_opts'],
        help='Options passed to the Pisces program.')

    parser.add_argument('--genus',
        default=def_configs['ref_genome_genus'],
        help='Genus of the reference genome.')
    parser.add_argument('--species',
        default=def_configs['ref_genome_species'],
        help='Species of the reference genome.')
    parser.add_argument('--build',
        default=def_configs['ref_genome_build'],
        help='Build of the reference genome.')

    args = parser.parse_args()

    sys.exit(main(args))

