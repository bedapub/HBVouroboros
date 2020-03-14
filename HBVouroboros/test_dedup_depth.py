import unittest

import depth

infile = '../testdata/test-dup-depth.txt'
outfile = '../testdata/test-dup-depth-out.txt'

depth.dedup_file(infile, outfile)
