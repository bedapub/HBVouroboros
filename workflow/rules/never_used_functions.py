#TODO: Verify they are used nowhere
def split_FASTA(infile, outdir=None, prefix=''):
    """	Split sequences in a FASTA file into separate files
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
        outfile = join(outdir, prefix + seq.id.replace('|', '_') + '.fasta')
        SeqIO.write(seq, outfile, 'fasta')
        count += 1
    
    return count

def get_simplified_id(desc):
    """ 
    Create a new id for reference genome that contains genotype and accession
   
    Args:
       desc (str): The input description
    Returns:
       str: a new id
    """

    strain = desc.split(" ")[0].split("|")[2]
    acc = strain.split("_")[0]
    gt = strain.split("-")[1]
    id = "{}|{}".format(gt, acc)
    
    return id

def set_samp_anno(perform_sim):
	"""Sets the appropriate smaples annotation file
	based on whether the pipeline is to be run with 
	simulated data, as specified by the 'do_sim' config 
	parameter"""
    
    #TODO: Take a look at the functionality
	if config['doSim'] == True:
                
		if perform_sim == True:
			sample_annotation = config['sample_annotation_sm']
			genomeId = config['genomeId']
			sampNum = config['sampNum']
			pairedEndDist = config['pairedEndDist']
			readLen = config['readLen']
			mt = config['mt']
			mp = config['mp']
			mtPos = config['mtPos']

			stream = os.system ('python RNAsim2/bin/RNAsim.py ' + ' '+" '"+genomeId+"' "+ "' "+sampNum+"' "+" '"+ pairedEndDist+"' "+" '"+ readLen+"' "+" '"+ mt+"' "+" '"+ mp+"' "+" '"+ mtPos+"' ")
                        
		return config['sample_annotation_sm']
	else:
		return config['sample_annotation']
