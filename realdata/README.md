Test mapping of real data using duplicated genomes
===

```bash
for id in 2 3 4; do
  bowtie2 -p 2  --no-mixed --no-discordant --sensitive -x ../HBV_refgenomes/HBV_refgenomes_dup_BOWTIE2/HBV_refgenomes_dup_BOWTIE2 2>log -1 HBV-DMSO-6d-"$id".unmapped_mate1.gz -2 HBV-DMSO-6d-"$id".unmapped_mate2.gz | samtools view -Sbh - > HBV-DMSO-6d-"$id".bam
done
```

Collect fastq files

```bash
../../fastqCollector/collect-unmapped-pe-fastq.bash . > sampleAnnotation.txt
```
