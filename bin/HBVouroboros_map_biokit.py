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
    indir = args.biokit_dir
    refgenomes_dir = args.refgenomes_dir

    outdir = args.outdir
    if outdir is None:
        outdir = os.path.join(indir, 'HBVouroboros')
    if os.path.exists(outdir):
        makedirs(outdir, mode=0x775, exist_ok=True)

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

    status = snakemake.snakemake(
        snakefile=align_snakefile,
        cluster='sbatch -t {cluster.time} -c {cluster.cpu} '
                  '-N {cluster.nodes} --mem={cluster.mem} '
                  '--ntasks-per-node={cluster.ntasks_per_node}',
        cluster_config=align_clusterfile,
        cores=64, nodes=64, local_cores=4,
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
       help = 'The directory of reference genomes')
    parser.add_argument('biokit_dir',
        help = 'An output directory of biokit')
    parser.add_argument('--outdir', '-o',
        help='output directory: if missing, a folder "HBVouroboros" '
             ' will be created in the biokit output directory')

    args = parser.parse_args()

    sys.exit(main(args))

