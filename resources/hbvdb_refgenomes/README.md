Reference genomes with features
===
Jitao David Zhang, June 19th, 2024.

I used the following scripts to download reference genomes from HBVdb on June 19th, 2024, and used the `seqret` tool from the EMBOSS suite to convert EMBL format files into GenBank files.

```bash
for strain in `cut -f 2 hbvdv_reference.tsv | awk 'NR>1'`; do wget https://hbvdb.lyon.inserm.fr/tmp/hbvdb_dat/${strain}/${strain}_entry.txt; done
ml load EMBOSS
for file in *_entry.txt; do genotype=`grep "genotype" $file | head -n 1 | sed 's/.*genotype //g' | sed 's/^\([A-J]\).*/\1/g'`; seqret ${file} -feature -osformat genbank -osname Genotype${genotype}_${file/_entry.txt/} -osextension gb -auto; done
```
