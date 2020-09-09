import unittest
from HBVouroboros import depth
import pandas as pd
import filecmp

class TestDedup(unittest.TestCase):

    def test_dedup(self):
      df = pd.DataFrame({
              '#CHR': ['chr1'] * 8,
              'POS' : range(1, 9),
              'S1' : [1, 4, 5, 6, 9, 1, 3, 5],
              'S2' : [3, 4, 5, 4, 3, 6, 3, 2], 
              'S3' : [5, 6, 3, 3, 4, 3, 2, 7]
              })
      df_dedup = depth.dedup(df)
      exp_dedup = pd.DataFrame({
          '#CHR': ['chr1'] * 4,
          'POS' : range(1,5),
          'S1' : [10, 5, 8, 11],
          'S2' : [6, 10, 8, 6],
          'S3' : [9, 9, 5, 10]
          })
      self.assertTrue(df_dedup.equals(exp_dedup))

    def test_dedup_file(self):
      infile = '../testdata/depth/test-dup-depth.txt'
      outfile = '../testdata/depth/test-dup-depth-out.txt'
      expfile = '../testdata/depth/test-dup-depth-exp-out.txt'

      depth.dedup_file(infile, outfile)
      self.assertTrue(filecmp.cmp(outfile, expfile))

if __name__ == '__main__':
    unittest.main()
