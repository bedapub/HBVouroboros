#!/usr/bin/env python3

import os.path

def collect_gene_coverage(coverage_files, outfile):
    """Collect gene coverage files into a GCT outfile

    Args:
        coverage_files (str): a list of coverage files exported by 
            'bedtools coverage'. They must be generated from the same
            GFF3 files.
        outfile (str): Output file name
    Returns:
        int: files collected
    """

    nsample = len(coverage_files)
    outf = open(outfile, 'w')

    ## coverage file name pattern (the logic is fragile - to be factored)
    ## infref_genome_{sample}_gene_coverage.tsv
    sample_names = [os.path.basename(f).
                      replace('infref_genome_', '').
                      replace('_feature_coverage.tsv', '')
                      for f in coverage_files]

    cov0 = coverage_files[0]
    genes = []
    descs = []

    with open(cov0, 'r') as f:
        for line in f:
            if line[0]=='#':
                next
            fl = line.rstrip().split('\t')
            if fl[2]=='gene':
                chrom = fl[0]
                gene = fl[8].replace("gene=","")
                genes.append(gene)
                descs.append(chrom)

    ngene = len(genes)

    outf.write('#1.2\n');
    outf.write('{}\t{}\n'.format(ngene, nsample))
    outf.write('Name\tDescription\t{}\n'.format(
                '\t'.join(sample_names)))

    counts=[[0]*ngene]*nsample

    for ind, cov in enumerate(coverage_files):
        cov_counts = []
        with open(cov, 'r') as f:
            for line in f:
                if line[0]=='#':
                    next
                fl = line.rstrip().split('\t')
                if fl[2]=='gene':
                    cov_counts.append(fl[9])
 
        counts[ind] = cov_counts

    for i in range(ngene):
        gene_counts = [vec[i] for vec in counts]
        curr_line = genes[i] + '\t' + \
                    descs[i] + '\t' + \
                    '\t'.join(gene_counts) + '\n'
        outf.write(curr_line)

    outf.close()
