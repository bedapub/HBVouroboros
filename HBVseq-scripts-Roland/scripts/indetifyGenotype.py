f = open(snakemake.input[0], "r")
l = f.read().split("\n")
f.close()
GenomeList = sorted(set([x.split("\t")[1][3:] for x in l[:-1] if (x[0:3] == "@SQ")]))
GenomeListMappedReads = [x.split("\t")[2] for x in l[:-1] if (x[0] != "@" and x.split("\t")[2] not in ["*", "0"])]
GenomeCount = [GenomeListMappedReads.count(x) for x in GenomeList]

gcd = (dict(zip(sorted(set(x for x in GenomeList if (GenomeListMappedReads.count(x)>0))), [GenomeListMappedReads.count(x) for x in GenomeList if (GenomeListMappedReads.count(x)>0)])))
GTL = sorted(set(x.split("|")[1] for x in GenomeList))
GTLMP = [x.split("|")[1] for x in GenomeListMappedReads]
GTC = [GTLMP.count(x) for x in GTL]
gtcd = (dict(zip(GTL, GTC)))
ou = open(snakemake.output[0], "w")
ou.write("Genotype:\tGenome:\tSequence:\n"+str(sorted(gtcd.items(), key=operator.itemgetter(1))[-1][0])+"\t"+
str(sorted(gcd.items(), key=operator.itemgetter(1))[-1][0])+l[l.index(">"+sorted(gcd.items(),
key=operator.itemgetter(1))[-1][0])+1]+"\n"+str(sorted(gtcd.items(), key=operator.itemgetter(1))[-1][1])+"\t"+
str(sorted(gcd.items(), key=operator.itemgetter(1))[-1][1]))
ou.close()
