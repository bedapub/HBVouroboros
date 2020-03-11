#!/usr/bin/env python3

import download_HBV_refgenomes
import dup_and_conc_fasta

download_HBV_refgenomes.download_HBV_refgenomes("test.fasta")
dup_and_conc_fasta.dup_and_conc_fasta("test.fasta", "test-dupconc.fasta")

