#!/usr/bin/env python3

import argparse
import sys
from os import makedirs
import os.path
import snakemake
import pkg_resources


align_snakefile = pkg_resources.resource_filename('HBVouroboros', 
    'align_reads/Snakefile')
align_clusterfile = pkg_resources.resource_filename('HBVouroboros',  
    'config/cluster.json')

if not os.path.exists(align_snakefile):
    raise Exception('align_snakefile not found: ' + align_snakefile)
if not os.path.exists(align_clusterfile):
    raise Exception('align_clusterfile not found:' + align_clusterfile)

def main(args):
    refgenomes_dir = os.path.realpath(args.refgenomes_dir)
    sample_file = os.path.realpath(args.sample_annotation_file)
    local = args.local
   
    if not local:
        cluster_comm = ('sbatch -t {cluster.time} -c {cluster.cpu} ' 
                       '-N {cluster.nodes} --mem={cluster.mem} ' 
                       '--ntasks-per-node={cluster.ntasks_per_node}')
        cluster_config = align_clusterfile
    else:
        cluster_comm = None
        cluster_config = None

    outdir = args.outdir
    if os.path.exists(outdir):
        makedirs(outdir, mode=0o775, exist_ok=True)

    status = snakemake.snakemake(
        snakefile=align_snakefile,
        cluster=cluster_comm,
        cluster_config=cluster_config,
        cores=128, nodes=128, local_cores=4,
        config={
            'sample_annotation': sample_file,
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
        description = 'Run HBVouroboros using an sample annotation file as input')

    parser.add_argument('refgenomes_dir',
       help = 'The directory of reference genomes')
    parser.add_argument('sample_annotation_file',
        help = 'A sample annotation file, a tab-delimited file with a header line. '
               'The first column contains sample IDs, the third column forward read files, and '
               'the fourth column reverse read files. Currently, single-read files are not tested.')
    parser.add_argument('--outdir', '-o',
        default='./HBVouroboros',
        help='output directory: if missing, a folder "HBVouroboros" '
             ' will be created in the current directory')
    parser.add_argument('-l', '--local', 
        action='store_true',
        help='If given, the commands will be executed locally.'
             'Otherwise, they will be submitted to the cluster with sbatch')

    args = parser.parse_args()

    sys.exit(main(args))

