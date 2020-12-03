import pytest
import RNAsim
import argparse
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


pytest.testSeqObject = SeqRecord(
    Seq("tttcacagctttccaacaagccctacaagatcccagagt"),
    id="gnl|hbvnuc|GQ924620_FT00000_C-C",
    name="HokC",
    description="gnl|hbvnuc|GQ924620_FT00000_C-C Feature FT:source AB076679_FT00000_P-A of [Viruses] Hepatitis B Virus genotype A (isolate HBV-Mala35). test fragment. [Hepadnaviridae] (length=39 residues).")

pytest.testSeqObject_dup = SeqRecord(
    Seq("tttcacagctttccaacaagccctacaagatcccagagttttcacagctttccaacaagccctacaagatcccagagt"),
    id="gnl|hbvnuc|GQ924620_FT00000_C-C|DupConc",
    name="HokC",
    description="gnl|hbvnuc|GQ924620_FT00000_C-C Feature FT:source AB076679_FT00000_P-A of [Viruses] Hepatitis B Virus genotype A (isolate HBV-Mala35). test fragment. [Hepadnaviridae] (original length=39 residues). Duplicated and concatenated (final length:78)")

@pytest.mark.parametrize("test_main_input,expected", [
    #correct required and optional arguments
    (['gnl|hbvnuc|GQ924620_FT00000_C-C', '2' ,'10' ,'3', '--fixedpos=56 66'], True),
    #readlength larger than fragment length
    (['gnl|hbvnuc|GQ924620_FT00000_C-C', '1', '4', '10'], False),
    #Wrong genotype name
    (['gnl|hbvnuc|AB076679_FT00000', '6' ,'10' ,'3'], False),
    #the number of specified positios not matching the number of paired reads required
    (['gnl|hbvnuc|GQ924620_FT00000_C-C', '2' ,'10' ,'3', '--fixedpos=56'], False)])

def test_main(test_main_input, expected):
    assert RNAsim.main(test_main_input) == expected

copy_testSeqObject = pytest.testSeqObject
def test_dup_and_conc():
    outcome = RNAsim.dup_and_conc(copy_testSeqObject)
    #Check duplication
    assert str(outcome.seq) == str(pytest.testSeqObject_dup.seq)
    #Check concatanation
    assert str(outcome.seq)[int(len(outcome.seq)/2)] == str(pytest.testSeqObject_dup.seq)[int(len(outcome.seq)/2)] == str(pytest.testSeqObject_dup.seq)[0]
    #Check id updating
    assert str(outcome.id) == str(pytest.testSeqObject_dup.id)
    #Check description updating
    assert outcome.description == pytest.testSeqObject_dup.description

def test_create_readPairs():
    #Using the manual specification of left hand starting point, verify the three (left, right, full fragment) generated sequencese
    leftSequences = []
    rightSequences = []
    fullFragment = []
    pairedEndDist = 8
    readLn = 4
    leftEnd = 12
    RNAsim.create_readPairs(pairedEndDist, readLn, leftEnd-1, pytest.testSeqObject_dup, leftSequences, rightSequences, fullFragment,1)
    assert str(rightSequences[0].seq) == "TTGT"
    assert str(leftSequences[0].seq) == "TCCA"
    assert str(fullFragment[0].seq) == "TCCAACAA"

    #test wether the reads span across the orginal genome en
    leftSequences = []
    rightSequences = []
    fullFragment = []
    leftEnd = 37
    RNAsim.create_readPairs(pairedEndDist, readLn, leftEnd-1, pytest.testSeqObject_dup, leftSequences, rightSequences, fullFragment,1)
    assert str(rightSequences[0].seq) == "TGAA"
    assert str(leftSequences[0].seq) == "AGTT"
    assert str(fullFragment[0].seq) == "AGTTTTCA"

copy_seq = str(pytest.testSeqObject.seq)
def test_point_mutate():
    #RNAsim.point_mutate(pytest.testSeqObject,15)
    outF = []
    makeSummary = False
    assert RNAsim.point_mutate(copy_testSeqObject, 15, outF, makeSummary) == True
    assert copy_seq[15] != str(copy_testSeqObject.seq[15])
    assert copy_seq[14] == str(copy_testSeqObject.seq[14])
    assert copy_seq[16] == str(copy_testSeqObject.seq[16])
