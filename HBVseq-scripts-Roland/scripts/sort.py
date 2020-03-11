import gzip
ma1 = gzip.open(snakemake.input[0], "r")
ma2 = gzip.open(snakemake.input[1], "r")
def takeFE(x):
    return(x[0])
l1 = []
for line in ma1:
    l1.append([str(line)[2:-7],str(ma1.readline())[2:-3],str(ma1.readline())[2:-3],str(ma1.readline())[2:-3]])
l2 = []
for line in ma2:
    l2.append([str(line)[2:-7],str(ma2.readline())[2:-3],str(ma2.readline())[2:-3],str(ma2.readline())[2:-3]])
ma1.close()
ma2.close()

l1.sort(key = takeFE)
l2.sort(key = takeFE)
ou1 = open(snakemake.output[0], "w")
ou2 = open(snakemake.output[1], "w")
a = snakemake.params[0]
for i in l1[:a]:
    ou1.write(i[0]+"\t"+"00"+"\n")
    ou1.write(i[1]+"\n")
    ou1.write(i[2]+"\n")
    ou1.write(i[3]+"\n")

for i in l2[:a]:
    ou2.write(i[0]+"\t"+"00"+"\n")
    ou2.write(i[1]+"\n")
    ou2.write(i[2]+"\n")
    ou2.write(i[3]+"\n")

ou1.close()
ou2.close()
