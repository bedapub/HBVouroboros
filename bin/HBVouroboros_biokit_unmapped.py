#!/usr/bin/env python3

import argparse
import sys
import os.path
import snakemake
import pkg_resources

align_snakefile = pkg_resources.resource_filename('HBVouroboros', 
    'align_reads/Snakefile')
align_clusterfile = pkg_resources.resource_filename('HBVouroboros',  
    'align_reads/config/cluster.yaml')

def main(args):
    indir = args.biokit_dir
    outdir = args.outdir
    if outdir == '':
       outdir = os.path.join(indir, 'HBVouroboros')
    sample_annotation = os.path.join(indir, 'samples.txt')

    status = snakemake.snakemake(
        snakefile = align_snakefile,
        cluster = 'sbatch -t {cluster.time} -c {cluster.cpu} '
                  '-N {cluster.nodes} --mem={cluster.mem} '
                  '--ntasks-per-node={cluster.ntasks_per_node}',
        cluster_config=align_clusterfile,
        cores=64, nodes=64, local_cores=4,
        config={'sample_annotation': sample_annotation},
        workdir=outdir,
        restart_times=3,
        printshellcmds=True)

    if status: # translate "success" into shell exit code of 0
        return 0
    return 1

 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Run HBVouroboros using a biokit outdir')

    parser.add_argument('biokit_dir',
        help = 'An output directory of biokit')
    parser.add_argument('--outdir', '-o',
        help='output directory: if missing, a folder "HBVouroboros" '
             ' will be created in the biokit output directory')

    args = parser.parse_args()

    sys.exit(main(args))

