#!/usr/bin/env python3

import os.path


def gene_coverage(coverage_file, outfile):
    """Report gene coverage

    Args:
        coverage_file (str): output file of 'bedtools coverage'
    Returns:
        int: genes written
    """
    
    ngenes = 0;
    out = open(outfile, 'w')
    out.write('Gene\tChrom\tcount\n')
    with open(coverage_file, 'r') as f:
        for line in f:
            fl = line.split('\t')
            if(fl[2]=='gene'):
                chrom = fl[0]
                gene = fl[8].replace("gene=","")
                count = fl[9]
                out.write('{}\t{}\t{}\n'.format(gene, chrom, count))
                ngenes += 1

    out.close()
    return(ngenes)

def collect_gene_coverage(coverage_files, outfile):
    """Collect gene coverage files into a GCT outfile

    Args:
        coverage_files (str): a list of coverage files
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
                      replace('_gene_coverage.tsv', '')
                      for f in coverage_files]

    cov0 = coverage_files[0]
    with open(cov0, 'r') as f:
        lines = f.readlines()
        ngene = len(lines)-1
        outf.write('#1.2\n');
        outf.write('{}\t{}\n'.format(ngene, nsample))
        outf.write('Name\tDescription\t{}\n'.format(
                '\t'.join(sample_names)))

    outf.close()
