#pat = "/home/roland/HBV/data/one_time_genomes/HBV-genotyped-completeGenomes-20190306.fasta"
unduplicated = open(snakemake.input[0], "r")
#pat2 = "/home/roland/HBV/data/duplicated_genomes/HBV-duplicatedGenomes.fasta"
duplicated = open(snakemake.output[0], "w")
lengths = open(snakemake.output[1], "w")
l = []
for line in unduplicated:
    tmp = unduplicated.readline()[:-1]
    if (line not in l):
        duplicated.write(line)

        duplicated.write(2*tmp+"\n")
        lengths.write(line+"\t"+str(len(tmp))+"\n")
        l.append(line)
lengths.close()
unduplicated.close()
duplicated.close()
