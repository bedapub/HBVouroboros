import requests
import logging
import os
from workflow.rules.python.common import dup_and_conc_FASTA

def download_refgenomes():
    refgenomes_response = requests.get("https://hbvdb.lyon.inserm.fr/data/references/hbvdbr.fas")
    open("resources/ref/HBV_refgenomes.fasta", "wb").write(refgenomes_response.content)

    dup_and_conc_FASTA("resources/ref/HBV_refgenomes.fasta", "resources/ref/HBV_refgenomes_dup.fasta")
    print("Generating bowtie2 index of duplicated HBV reference genomes")
    stream = os.popen("mkdir resources/ref/HBV_refgenomes_dup_BOWTIE2; bowtie2-build resources/ref/HBV_refgenomes_dup.fasta resources/ref/HBV_refgenomes_dup_BOWTIE2; touch resources/ref/HBV_refgenomes_dup_BOWTIE2")
    logging.info(stream.read())

def download_all_genomes():
    refgenomes_response = requests.get("https://hbvdb.lyon.inserm.fr/data/nucleic/fasta/all_Genomes.fas")
    open("resources/ref/HBV_allgenomes.fasta", "wb").write(refgenomes_response.content)

    stream = os.popen("makeblastdb -in resources/ref/HBV_allgenomes.fasta -title \"HBVdb genomes\" -dbtype nucl")
    logging.info(stream.read())

def download_files():
    logging.basicConfig(filename="logs/ref/build_refgenomes_python.log", level=logging.INFO)
    download_refgenomes()
    download_all_genomes()

download_files()