Test mapping using duplicated genomes
===

```bash
for id in 1 2 3 4; do
  bowtie2 -p 2  --no-mixed --no-discordant --sensitive -x ../HBV_refgenomes/HBV_refgenomes_dup_BOWTIE2/HBV_refgenomes_dup_BOWTIE2 2>log -1 testfile-"$id"_1.fq.gz -2 testfile-"$id"_2.fq.gz | samtools view -Sbh - > testfile-"$id".bam
done
```

```bash
for id in 2 3 4; do
  bowtie2 -p 2  --no-mixed --no-discordant --sensitive -x ../HBV_refgenomes/HBV_refgenomes_dup_BOWTIE2/HBV_refgenomes_dup_BOWTIE2 2>log -1 HBV-DMSO-6d-"$id".unmapped_mate1.gz -2 HBV-DMSO-6d-"$id".unmapped_mate2.gz | samtools view -Sbh - > HBV-DMSO-6d-"$id".bam
done

```
