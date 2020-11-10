#!/usr/bin/env python3

from os.path import join, exists, isdir, isfile, expanduser, \
    realpath, basename, dirname
from pandas import read_table
from Bio import SeqIO
from Bio import Entrez
from BCBio import GFF


def biokit_sample_annotation_filename(biokit_outdir):
    return(join(biokit_outdir, 'samples.txt'))


def parse_sample_annotation(sample_annotation_file):
    """ Parse biokit sample annotation file into sample names and FASTQ dicts

    Args:
        sample_annotation_file (str): A Biokit sample annotation file
    Returns:
        A tuple of three elements: samples, fq1dict (dict of FASTQ1
        files indexed by sample names), fq2dict (dict of FASTQ2 files
        indexed by sample names).
    """

    annotation = read_table(sample_annotation_file)
    samples = annotation.iloc[:, 0]
    fq1s = annotation.iloc[:, 2]
    fq2s = annotation.iloc[:, 3]

    fq1dict = dict(zip(samples, fq1s))
    fq2dict = dict(zip(samples, fq2s))

    return(samples, fq1dict, fq2dict)


def biokit_unmapped_sample_annotation(biokit_outdir, outfile):
    """ Get sample annotation from a biokit output directory

    Args:
        biokit_outdir (str): An output directory of the biokit pipeline
        outfile (str): Output file name of sample annotation
    Returns:
        number of samples
    """

    biokit_outdir = realpath(expanduser(biokit_outdir))
    if not isdir(biokit_outdir):
        raise Exception(
            'Directory {} does not exist'.format(biokit_outdir))
    infile = biokit_sample_annotation_filename(biokit_outdir)
    if not exists(infile):
        raise Exception(
           'sample annotaiton file ({}) not found'.format(infile))

    fout = open(outfile, 'w')
    with open(infile, 'r') as fin:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            lsplit = line.rstrip().split('\t')
            fastq_format = join(
              biokit_outdir,
              'unmapped',
              '.'.join([
                   lsplit[0],
                   'unmapped_mate{}',
                   'gz']))
            f1 = fastq_format.format(1)
            f2 = fastq_format.format(2)
            if not isfile(f1):
                raise Exception('File {} does not exist'.format(f1))
            if not isfile(f2):
                raise Exception('File {} does not exist'.format(f2))
            lsplit[2] = f1
            lsplit[3] = f2
            newline = '\t'.join(lsplit)
            fout.write(newline)
            fout.write('\n')

    fout.close()


def collect_gene_coverage(coverage_files, outfile, feat_type='gene'):
    """Collect gene coverage files into a GCT outfile

    Args:
        coverage_files (str): a list of coverage files exported by
            'bedtools coverage'. They must be generated from the same
            GFF3 files.
        outfile (str): Output file name
        feat_type (str): feature type, 'gene' or 'CDS'
    Returns:
        int: files collected
    """

    nsample = len(coverage_files)
    outf = open(outfile, 'w')

    # coverage file name pattern (the logic is fragile - to be factored)
    # infref_genome_{sample}_gene_coverage.tsv
    sample_names = [basename(f).
                    replace('infref_genome_', '').
                    replace('_feature_coverage.tsv', '')
                    for f in coverage_files]

    cov0 = coverage_files[0]
    genes = []
    descs = []

    with open(cov0, 'r') as f:
        for line in f:
            if line[0] == '#':
                next
            fl = line.rstrip().split('\t')
            if fl[2] == feat_type:
                chrom = fl[0]
                if feat_type == 'gene':
                    gene = fl[8].replace("gene=", "")
                elif feat_type == 'CDS':
                    # this assumes that the ID is in the second
                    # which is not always the case - to be fixed
                    gene = fl[8].split(';')[1].split('=')[1]
                else:
                    raise Exception('not supported feat_type: CDS/gene only')
                genes.append(gene)
                descs.append(chrom)

    ngene = len(genes)

    outf.write('#1.2\n')
    outf.write('{}\t{}\n'.format(ngene, nsample))
    outf.write('Name\tDescription\t{}\n'.format(
                '\t'.join(sample_names)))

    counts = [[0]*ngene]*nsample

    for ind, cov in enumerate(coverage_files):
        cov_counts = []
        with open(cov, 'r') as f:
            for line in f:
                if line[0] == '#':
                    next
                fl = line.rstrip().split('\t')
                if fl[2] == feat_type:
                    cov_counts.append(fl[9])

        counts[ind] = cov_counts

    for i in range(ngene):
        gene_counts = [vec[i] for vec in counts]
        curr_line = genes[i] + '\t' + \
            descs[i] + '\t' + \
            '\t'.join(gene_counts) + '\n'
        outf.write(curr_line)

    outf.close()


def dedup(depth):
    """ Deduplicate depth output of samtools

    Args:
        depth (pandas.DataFrame): first two columns are chromsomes and
            depths, respectively, and the rest columns are samples.
            Positions are duplicated
    Returns:
        pandas.DataFrame: a new DataFrame containing position aggregated
            values
    """

    pos = depth.iloc[:, 1]
    duplen = len(pos)

    if duplen % 2 != 0:
        raise Exception("the input file has odd-number positions")
    olen = int(duplen/2)

    out = depth.iloc[:olen, :].copy()
    for i in range(olen, duplen):
        out.iloc[i-olen, 2:] = out.iloc[i-olen, 2:] + depth.iloc[i, 2:]

    # bam file name pattern (the logic is fragile - to be factored)
    # infref_bam/Sample4.sorted.bam
    out.rename(mapper=lambda colname:
               basename(colname).replace('.sorted.bam', ''),
               axis='columns', inplace=True)

    return(out)


def dedup_file(infile, outfile):
    """Perform dedup on the output file of samtools depth
    Args:
        infile (str): depth file exported by samtools *with the header*
        outfile (str): outfile file with deduplicated depths

    Returns:
        value of pandas.DataFrame.to_csv
    """

    depth = read_table(infile)
    dedup_res = dedup(depth)
    res = dedup_res.to_csv(outfile, sep='\t', index=False)
    return(res)


def get_infref_acc(blast_tab_file):
    """Get accession of the inferred reference strain

    Args:
        blast_tab_file (str): Tabular output file of BLAST (outfmt=6)

    Returns:
        str: accession number
    """

    with open(blast_tab_file, 'r') as f:
        return(f.read().split('\t')[1])


def get_infref_gb_acc(blast_tab_file):
    """Get GenBank accession of the inferred reference strain

    Args:
        blast_tab_file (str): Tabular output file of BLAST (outfmt=6)

    Returns:
        str: genbank accession number
    """

    acc = get_infref_acc(blast_tab_file)
    res = acc.split("|")[2].split("_")[0]
    return(res)


def download_gb(acc, outfile):
    """Download GenBank file with the given accession number

    Args:
        acc (str): GenBank acession number, example:'KC774468'
        outfile (str): output GenBank file

    Returns:
        int: number of records written as an integer.
    """

    res = -1
    Entrez.email = "jitao_david.zhang@roche.com"
    with Entrez.efetch(
      db="nuccore", rettype="gb", retmode="text", id=acc
    ) as handle:
        seq_record = SeqIO.read(handle, "gb")

    res = SeqIO.write(seq_record, outfile, 'gb')
    return(res)


def write_seq_by_acc(infile, acc, outfile):
    """Fetch sequence by accession number and write it to FASTA file

    Args:
        infile (str): FASTA file of genomes
        acc (str): Accession number of the inferred reference strain
        outfile (str): output FASTA file
    Returns:
        bool: if found or not
    """

    found = False
    seqs = SeqIO.parse(infile, 'fasta')
    for seq in seqs:
        if seq.id == acc:
            found = True
            SeqIO.write(seq, outfile, 'fasta')

    return(found)


def gb2gff(infile, outfile):
    """Translate GenBank file to GFF3 file. TODO: the procedure now does not
    handle join correctly

    Args:
        infile (str): input GenBank file
        outfile (str): output GFF3 file
    Returns:
        Number of records written
    """

    gb_handle = open(infile, 'r')
    gff_handle = open(outfile, 'w')
    res = GFF.write(SeqIO.parse(gb_handle, "gb"), gff_handle)
    gff_handle.close()
    return(res)


def sort_FASTA_by_length(infile, outfile):
    """Sort FASTA sequence by length descendingly

    Args:
        infile (str): FASTA file
        outfile (str): Output FASTA file
    return:
        None
    """

    records = list(SeqIO.parse(infile, "fasta"))
    records.sort(key=lambda r: -len(r))
    SeqIO.write(records, outfile, "fasta")


def first_accession(fastafile):
    """Get accession of the first record in FASTA file

    Args:
        fastafile (str): FASTA file
    Returns:
        str : accession number of the first record
    """

    res = ''
    seqs = SeqIO.parse(fastafile, 'fasta')
    for seq in seqs:
        res = seq.id
        break

    return(res)


def dup_gff(dup_fasta, ingff, outgff):
    """Make GFF files for duplicated genome

    Args:
        dup_fasta (str): duplicated FASTA
        ingff (str): input GFF file
        outgff (str): output GFF file
    Returns:
        None
    """

    acc = first_accession(dup_fasta)
    out_handle = open(outgff, 'w')
    with open(ingff, 'r') as f:
        for line in f:
            if line[0] == '#':
                out_handle.write(line)
            else:
                fl = line.split('\t')
                fl[0] = acc
                out_handle.write('\t'.join(fl))
            out_handle.write('\n')

    out_handle.close()
    return(None)


def get_simplified_id(desc):
    """ Make a new id for reference genome that contains genotype and accession

    Args:
       desc (str): The input description

    Returns:
       str: a new id
    """

    strain = desc.split(" ")[0].split("|")[2]
    acc = strain.split("_")[0]
    gt = strain.split("-")[1]
    id = "{}|{}".format(gt, acc)
    return(id)


def dup_and_conc(record):
    """Duplicate the sequence and concatenate the original and duplicated
    sequence, and append a text label to the id and the description

    Args:
        record (Bio.SeqRecord): A SeqRecord object
    Returns:
        Bio.SeqRecord
    """

    record.seq = record.seq*2
    old_desc = record.description
    new_desc = old_desc.replace("length=", "original length=")
    new_desc += ' Duplicated and concatenated \
            (final length:{})'.format(len(record.seq))
    record.description = new_desc
    record.id += '|DupConc'
    return(record)


def dup_and_conc_FASTA(infile, outfile):
    """Duplicate sequences in a FASTA file, concatenate the original
       with the duplicates, and write the concatenated sequences

    Args:
        infile (str): The input filename, pointing to a FASTA file
        outfile (str): The output filename, overwritten if the file exists.

    Returns:
        int: number of sequences
    """

    sequences = SeqIO.parse(infile, 'fasta')
    outseqs = (dup_and_conc(record) for record in sequences)
    res = SeqIO.write(outseqs, outfile, 'fasta')
    return(res)


def split_FASTA(infile, outdir=None, prefix=''):
    """Split sequences in a FASTA file into separate files
       output file name is given by the ids (with pipes replaced by underscore)

       Args:
           infile (str): The input FASTA file name
           outdir (str): The output directory. Default: input file folder
           prefix (str): Prefix to the output file name
       Returns:
           int: number of sequences
    """

    count = 0
    if outdir is None:
        outdir = dirname(infile)
    sequences = SeqIO.parse(infile, 'fasta')
    for seq in sequences:
        outfile = join(outdir,
                       prefix + seq.id.replace('|', '_') + '.fasta')
        SeqIO.write(seq, outfile, 'fasta')
        count += 1
