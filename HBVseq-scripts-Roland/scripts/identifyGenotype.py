import operator
import re
f = open(snakemake.input[0], "r")
l = f.read().split("\n")
#all genome lines
lll = [x for x in l[:-1] if x[0:3]=="@SQ"]
#all genomes and their length
gld = dict(zip([x.split("\t")[1][3:] for x in lll],[int(x.split("\t")[2][3:])/2 for x in lll]))
#all reads
ll = [x for x in l[:-1] if (x[0]!="@" and min(int(x.split("\t")[3]),
int(x.split("\t")[7]))<gld.get(x.split("\t")[2]))]
f.close()
#all genomes
GenomeList = sorted(set(x.split("\t")[1][3:] for x in lll))
#all genomes with reads
GenomeListMappedReads = [x.split("\t")[2] for x in ll if (x[0] != "@" and
 x.split("\t")[2] not in ["*", "0"])]
#count off reads per genome
GenomeCount = [GenomeListMappedReads.count(x) for x in GenomeList]
#all genomes with 95% off the most reads as keys and their reads as values
gcd = (dict(zip(sorted(
set(x for x in GenomeList if (GenomeListMappedReads.count(x)>=(0.95*max(GenomeCount))))),
[GenomeListMappedReads.count(x) for x in GenomeList if
(GenomeListMappedReads.count(x)>=(0.95*max(GenomeCount)))])))
#all Genotypes
GTL = sorted(set(x.split("|")[1] for x in GenomeList))

ou = open(snakemake.output[0], "w")
#orderd gcd in regards to reads
sgcd  =sorted(gcd.items(), key=operator.itemgetter(1), reverse = 1)


most_mapped_genomes_genotypes = [x[0].split("|")[1] for x in sgcd]
ou.write("Genotype:\n")
for gt in GTL:
    ou.write(gt+":\t"+str(most_mapped_genomes_genotypes.count(gt))+"\t")
ou.write("\nGenomes:\n")
for gt, readc in sgcd:
    ou.write(str(gt)+"\t"+str(readc)+"\n")

ou.close()
#
genomesf = open("data/duplicated_genomes/HBV.fasta","r")
#genome sequences
gtemp = {line[1:-1]:genomesf.readline() for line in genomesf}
genomesf.close()
#sample sequences
samples = set(x for x in ll if x.split("\t")[9] not in ["0","*"])

for genome in sgcd:
#all reads for a genome
    temp = [x.split("\t") for x in ll if x.split("\t")[2]==genome[0]]
    #which parts of the genome get read
    tmp =[range(int(x[3]), int(x[3])+int(x[8])) for x in temp if x[8][0]=="-"]
    #all base-positions which get read
    tmp2 =[i%int(gld.get(genome[0])) for x in temp for i in list(range(int(x[3]), int(x[3])+int(x[8])))


    for x in temp:
        #positioin in the genome
        i = int(x[3])
#position in the sample
        c = 0
        for match in re.findall("[0-9]+[X=IDMNSHPX]", x[5]):
            if (match[-1] in ["X","=","D"]):
                if (match[-1]=="X"):
                    print("Pointmutation(s) at: "+str(i%int(gld.get(genome[0])))+" till "+
                    str((i+int(match[:-1]))%int(gld.get(genome[0])))+" from: "+
                    gtemp.get(genome)[i:i+int(match[:-1])]+" to: "+[y.split("\t")[9] for y in samples if y.split("\t")[0] == x[0] and
                    y.split("\t")[8][0] == x[8][0]][0][c:c+int(match[:-1])])

                if (match[-1] == "D"):
                    print("Deletion(s) in the sample "+x[0]+" from: "+str(i%int(gld.get(genome[0])))+
                    " to : "+
                    str((i+int(match[:-1]))%int(gld.get(genome[0]))))
                i = i+int(match[:-1])
                c= c+int(match[:-1])
            if (match[-1] in ["I"]):
                print("Insertion(s) in the sample "+x[0]+" from: "+str(i%int(gld.get(genome[0])))+
                " to : "+
                str((i+int(match[:-1]))%int(gld.get(genome[0]))))
print(temp)
print(tmp)
ma1 = gzip.open("data/test/Z-2.unmapped_mate1.gz","r")
ma2 = gzip.open("data/test/Z-2.unmapped_mate2.gz","r")
for line in ma1:
