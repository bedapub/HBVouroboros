Test mapping using duplicated genomes
===

```bash
for id in 1 2 3 4; do
  bowtie2 -p 2  --no-mixed --no-discordant --sensitive -x ../HBV_refgenomes/HBV_refgenomes_dup_BOWTIE2/HBV_refgenomes_dup_BOWTIE2 2>log -1 testfile-"$id"_1.fq.gz -2 testfile-"$id"_2.fq.gz | samtools view -Sbh - > testfile-"$id".bam
done
```
