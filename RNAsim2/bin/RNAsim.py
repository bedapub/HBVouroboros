import Bio
import copy
import random
import sys
import re
import os
import os.path
from Bio import SeqIO, bgzf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import pathlib
from gzip import open as gzopen

def create_readPairs(pairedEndDist, readLen, leftEnd, record_dup, leftSequences, rightSequences, fullFragment,i):

    """Geneartes paired end reads, saves left and right read and the full fragment

    Args:
        pairedEndDist (int)
        readLen (int)
        noReads
        record_dup (Bio.SeqRecord): A SeqRecord object
        leftSequences (Bio.SeqRecord): A SeqRecord object
        rightSequences (Bio.SeqRecord): A SeqRecord object
        fullFragment (Bio.SeqRecord): A SeqRecord object
    Returns:
    none
    """

    rightEnd = leftEnd + pairedEndDist

    #the fragment from a radial genome
    if rightEnd > len(record_dup.seq):
        circSeq = record_dup.seq[leftEnd:len(record_dup.seq)] + record_dup.seq[0:rightEnd-len(record_dup.seq)]
        #print(circSeq)
    else:
        circSeq = record_dup.seq[leftEnd:rightEnd]
        #print(circSeq)

    #rightEnd = leftEnd + pairedEndDist

    leftHandread = SeqRecord(circSeq[0:readLen].upper(),id=str(i), description="")
    rightHandread = SeqRecord(circSeq[len(circSeq)-readLen:len(circSeq)].upper().reverse_complement(),id=str(i),description="")
    theFullFragment = SeqRecord(circSeq.upper(),id='F_'+str(i),description="")
    #asigns some arbitary id, labels etc for now
    #leftHandread = SeqRecord(record_dup.seq[leftEnd:leftEnd+readLen].upper(),id=str(i), description="")
    leftSequences.append(leftHandread)
    #rightHandread = SeqRecord(record_dup.seq[rightEnd-readLen:rightEnd].upper().reverse_complement(),id=str(i),description="")
    rightSequences.append(rightHandread)
    #theFullFragment = SeqRecord(record_dup.seq[leftEnd:rightEnd].upper(),id='F_'+str(i),description="")
    fullFragment.append(theFullFragment)




    #asigns some arbitary id, labels etc for now
    #leftHandread = SeqRecord(circSeq[0:readLen].upper(),id=str(i), description="")
    #leftSequences.append(leftHandread)
    #rightHandread = SeqRecord(circSeq[len(circSeq)-readLen:len(circSeq)].upper(),id=str(i),description="")
    #rightSequences.append(rightHandread)
    #theFullFragment = SeqRecord(circSeq.upper(),id='F_'+str(i),description="")
    #fullFragment.append(theFullFragment)

    #print(leftHandread.seq)
    #print(rightHandread.seq)


def point_mutate(record_copy,theMutationPos,outF, makeSummary):
    """Change the base in record_copy.seq at theMutationPos to another random base

    Args:
        record_copy (Bio.SeqRecord): A SeqRecord object
        theMutationPos: int
    Returns:
        Boolean
    """
    NAs = ["A", "T", "G", "C"]
    originalB = record_copy.seq[theMutationPos]
    mutated_seq = list(str(record_copy.seq))
    #identify the original base and replace in with another random base
    NAs = [ x for x in NAs if x!=originalB]
    mutated_seq[theMutationPos] = random.choice(NAs)
    #convert back the mutated sequqnce to str and replace the original sequence with it
    mutated_seq_final =""
    for x in mutated_seq:
        mutated_seq_final += x
    record_copy.seq = Seq(mutated_seq_final)
    if makeSummary == True:
        outF.write(str(originalB) + str(theMutationPos) + str(mutated_seq[theMutationPos]))
        outF.write("\n")

    return (True)

def dup_and_conc(record):
    """Duplicate the sequence and concatenate the original and duplicated
    sequence, and append a text label to the id and the description

    Args:
        record (Bio.SeqRecord): A SeqRecord object
    Returns:
        record (Bio.SeqRecord)
    """
    record.seq = record.seq*2
    old_desc = record.description
    new_desc = old_desc.replace("length=", "original length=")
    new_desc += ' Duplicated and concatenated (final length:{})'.format(len(record.seq))
    record.description = new_desc
    record.id += '|DupConc'
    return(record)

def main(args):

    parser = argparse.ArgumentParser(
        description='specify paramters to run a RNA paired-end sampling')
    parser.add_argument('genotypeId',
        help = 'the id of the genotype of choice')
        #e.g. gnl|hbvnuc|AB076679_FT00000_P-A
    parser.add_argument('noReads',
        help = 'number of paired end samples')
    parser.add_argument('pairedEndDist',
        help = 'the length of the full fragment')
    parser.add_argument('readLen',
        help = 'length of each of the paired reads')
    parser.add_argument('-fp', '--fixedpos',
        help="Specify the starting left-hand position of each read")
    parser.add_argument('-mt','--mutate', action="store_true",
        help="Introduce mutation to the reference genome")
    parser.add_argument('-rm', '--randmut',
        help="mutations will be preformed randomly")
    parser.add_argument('-mp', '--mutpos',
        help="specify the posisiton of each point mutation to tbe introduced")
    parser.add_argument('-p', '--percent',
        help="The percentrage of the samples that will be sampled from the genome after application fo specified mutations")
    args = parser.parse_args(args)



    #Find the path to the cloned RNAsim directory on the local machine
    thePath = str(pathlib.Path(__file__).parent.absolute())
    srcDir = thePath.replace('/bin', '')

    noReads = int(args.noReads)
    pairedEndDist = int(args.pairedEndDist)
    readLen = int(args.readLen)
    posArray = None
    mutPercent = 0
    if args.fixedpos:
        posArray = list(map(int, args.fixedpos.split(" ")))
        if len(posArray) != noReads:
            print("the number of specified positions does nto match the number of reads specified")
            return(False)
            sys.exit()

    if readLen > pairedEndDist:
        print("Bad input: read length is larger than fragment length")
        return(False)
        sys.exit()
    if os.path.isdir(srcDir) == False:
        print("The specified path does not exist")
        return (False)
        sys.exit()

    if args.mutate:
        if args.mutpos:
            mutPos = {}
            pos_list = args.mutpos.split()
            for i in range (0, len(pos_list)):
                pos_list[i]=int(pos_list[i])
            mutPos = pos_list
            if mutPos == None:
                print("You have chosen to manually specify mutation positions but have not provided a list of positions to be mutated")
                sys.exit()
        elif args.randmut:
            randMut= args.randmut
        else:
            print("You have chosen to simulate mutations but not provided the necessary argument to choose a random or deterministic mutation implementation")
            sys.exit()

    if args.percent:
        if args.mutate is None:
            print("mutation percentage specified but no mutation arguments are given")
            sys.exit()
        if args.percent >1 or arg.percent<0:
            print("The specified percentage must be a number between 0 and 1")
        mutpercent = args.percent

    #Find the input genome based on id  
    rootDir = srcDir
    
    oneSequence = SeqIO.parse(os.path.join(rootDir.replace('/RNAsim2','') + '/resources/ref/HBV_refgenomes.fasta'), 'fasta')
    record_copy = None
    for record in oneSequence:
        if record.id == args.genotypeId:

            record_copy = copy.deepcopy(record)
            orgRecord = copy.deepcopy(record)
            break
    if record_copy is None:
        print("No matching genome was found")
        return(False)
        sys.exit()
    else:
        if posArray != None:
            for element in posArray:
                if element > int(len(record_copy)):
                    print("Bad input: fragmenet length specified by fixedPos too small to accomiadate the fragment length")
                    return(False)
                    sys.exit()
    #Do genmoe mutation
    if args.mutate:
        outF = open(srcDir +"/output/mutation_summary.txt", "w")
        if args.randmut:
            #For the random mutation choose the given number of positions randomly
            for j in range(0, int(randMut)):
                point_mutate(record_copy, random.randint(0,len(record_copy.seq)) )

        if args.mutpos:
            for j in range(0, len(mutPos)):
                if int(mutPos[j]) > len(record_copy):
                    print("mutation position exceeds the size fo the genome")
                    sys.exit()
                makeSummary = True
                point_mutate(record_copy, mutPos[j], outF, makeSummary)
        outF.close()

    record_dup = dup_and_conc(record_copy)
    orgRecord_dup = dup_and_conc(orgRecord)


    leftSequences = []
    rightSequences = []
    fullFragment = []
    counter = 0

    for i in range(0,noReads):
        if args.percent:
            mutReads = int(noReads*mutpercent)
            counter = counter +1
            if counter <= mutReads:
                if posArray is None:
                    leftEnd = random.randint(0,len(orgRecord_dup.seq)/2)
                    create_readPairs(pairedEndDist, readLen, leftEnd, orgRecord, leftSequences, rightSequences, fullFragment, i)
                else:
                    create_readPairs(pairedEndDist, readLen, posArray[i]-1, orgRecord, leftSequences, rightSequences, fullFragment, i)
            else:
                if posArray is None:
                    leftEnd = random.randint(0,len(record_dup.seq)/2)
                    create_readPairs(pairedEndDist, readLen, leftEnd, record_copy, leftSequences, rightSequences, fullFragment, i)
                else:
                    create_readPairs(pairedEndDist, readLen, posArray[i]-1, record_copy, leftSequences, rightSequences, fullFragment, i)
        else:
            if posArray is None:
                leftEnd = random.randint(0,len(record_dup.seq)/2)
                create_readPairs(pairedEndDist, readLen, leftEnd, record_copy, leftSequences, rightSequences, fullFragment, i)
            else:
                create_readPairs(pairedEndDist, readLen, posArray[i]-1, record_copy, leftSequences, rightSequences, fullFragment, i)


        #Fake read quality
        leftSequences[i].letter_annotations["phred_quality"] = [50] * len(leftSequences[i].seq)
        rightSequences[i].letter_annotations["phred_quality"] = [50] * len(rightSequences[i].seq)
        fullFragment[i].letter_annotations["phred_quality"] = [50] * len(fullFragment[i].seq)

    with bgzf.BgzfWriter(os.path.join(srcDir + '/output/simSample-1_1.fastq.gz'), "wb") as outgz:
        SeqIO.write(sequences=leftSequences, handle=outgz, format="fastq")
    with bgzf.BgzfWriter(os.path.join(srcDir + '/output/simSample-1_2.fastq.gz'), "wb") as outgz:
        SeqIO.write(sequences=rightSequences, handle=outgz, format="fastq")

    sampAnnotation = open (os.path.join(srcDir+'/output/sampleAnnotation'), "w+")
    #sampAnnotation.truncate(0)
    sampAnnotation.write("#ID	GROUP	FASTQ1	FASTQ2\n")
    sampAnnotation.write("simSample	control	%s	%s" %(os.path.abspath(os.path.join(srcDir +'/output/simSample-1_1.fastq.gz')), os.path.abspath(os.path.join(srcDir +'/output/simSample-1_2.fastq.gz'))))
    sampAnnotation.close()


    return(True)


if __name__ == '__main__':
    main(sys.argv[1:])
    #sys.exit(main(args))
