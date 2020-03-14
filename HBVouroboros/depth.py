#!/usr/bin/env python3

import pandas as pd

def dedup(depth):
    """ Deduplicate depth output of samtools

    Args:
        depth (pandas.DataFrame): first two columns are chromsomes and 
            depths, respectively, and the rest columns are samples. 
            Positions are duplicated
    return:
        pandas.DataFrame: a new DataFrame containing position aggregated
            values
    """

    pos = depth.iloc[:, 1]
    duplen = len(pos)
    
    if duplen % 2!= 0:
        raise Exception("the input file has odd-number positions")
    olen = int(duplen/2)

    out = depth.iloc[:olen,:].copy()
    for i in range(olen, duplen):
        out.iloc[i-olen, 2:] = out.iloc[i-olen, 2:] + depth.iloc[i, 2:]

    return(out)

def dedup_file(infile, outfile):
    """Perform dedup on the output file of samtools depth
    Args:
        infile (str): depth file exported by samtools *with the header*
        outfile (str): outfile file with deduplicated depths
    return:
        value of pandas.DataFrame.to_csv
    """
    
    depth = pd.read_table(infile)
    dedup_res = dedup(depth)
    res = dedup_res.to_csv(outfile, sep='\t', index=False)
    return(res)
