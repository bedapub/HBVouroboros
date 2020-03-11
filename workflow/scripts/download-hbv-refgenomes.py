#!/usr/bin/env python3

import requests

def download_HBV_reference_genomes(filename="hbvdbr.fas"):
  """Download HBV reference genomes and save them into a file

  Args:
      filename (str): The output filename

  Returns:
      bool: whether the request was ok
  """

  url = 'https://hbvdb.lyon.inserm.fr/data/references/hbvdbr.fas'
  r = requests.get(url)
  rstr = r.content.decode("utf-8")
  with open(filename, 'w') as f:
      f.write(rstr)

  return(r.ok)
