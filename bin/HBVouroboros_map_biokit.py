#!/usr/bin/env python3

import argparse
import sys
from os import makedirs
import os.path
import snakemake
import pkg_resources

biokit_snakefile = pkg_resources.resource_filename('HBVouroboros',
    'biokit/Snakefile')
align_snakefile = pkg_resources.resource_filename('HBVouroboros', 
    'align_reads/Snakefile')
align_clusterfile = pkg_resources.resource_filename('HBVouroboros',  
    'align_reads/config/cluster.json')

if not os.path.exists(align_snakefile):
    raise Exception('align_snakefile not found')
if not os.path.exists(align_clusterfile):
    raise Exception('align_clusterfile not found')

def main(args):
    indir = os.path.realpath(args.biokit_dir)
    refgenomes_dir = os.path.realpath(args.refgenomes_dir)

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.join(indir, 'HBVouroboros')
    if not os.path.exists(outdir):
        makedirs(outdir, mode=0o775, exist_ok=True)

    unmapped_sample_annotation=os.path.join(outdir,
        'unmapped_samples.txt')
    biokit_status = snakemake.snakemake(
        snakefile=biokit_snakefile,
        config={
            'biokit_dir': indir,
            'output_file': unmapped_sample_annotation
        },
        workdir=outdir)
    if not biokit_status:
        raise Exception('Failed to derive sample annotation files '
            'for unmapped reads in directory {}'.format(indir))

    cluster_logs_dir = os.path.join(outdir, 'cluster-logs')
    makedirs(cluster_logs_dir, mode=0o775, exist_ok=True)
    cluster_out_pattern = os.path.join(cluster_logs_dir, 'slurm-%x-%j.out')
    cluster_err_pattern = os.path.join(cluster_logs_dir, 'slurm-%x-%j.err')

    if not args.local:
        cluster_comm = ('sbatch -t {cluster.time} -c {cluster.cpu} '
                       '-N {cluster.nodes} --mem={cluster.mem} '
                       '--ntasks-per-node={cluster.ntasks_per_node}'
                       ' -o ' + cluster_out_pattern + \
                       ' -e ' + cluster_err_pattern)
        cluster_config = align_clusterfile
    else:
        cluster_comm = None
        cluster_config = None

    status = snakemake.snakemake(align_snakefile,
        cluster=cluster_comm,
        cluster_config=cluster_config,
        cores=128, nodes=128, local_cores=4,
        config={
            'sample_annotation': unmapped_sample_annotation,
            'refgenomes_dir': refgenomes_dir
            },
        workdir=outdir,
        restart_times=3,
        printshellcmds=True)

    if status: # translate "success" into shell exit code of 0
        return 0
    return 1

 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Run HBVouroboros using a biokit outdir')

    parser.add_argument('refgenomes_dir',
       help = 'Path to the directory of reference genomes.')
    parser.add_argument('biokit_dir',
        help = 'Path to the output directory of the biokit pipeline')
    parser.add_argument('--outdir', '-o',
        help='output directory: if missing, a folder "HBVouroboros" '
             ' will be created in the biokit output directory')
    parser.add_argument('-l', '--local',
        action='store_true',
        help='If given, the commands will be executed locally.'
             'Otherwise, they will be submitted to the cluster with sbatch')

    args = parser.parse_args()

    sys.exit(main(args))

