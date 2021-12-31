rm(list = ls())
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

dd = as.data.frame(fread(args[1]))
opts = colnames(dd)
refgenome = strsplit(dd$`#CHROM`[1], split='\\|')[[1]][3]
png(file = "results/coverage/infref_genome_depth_mqc.png",
    width = 1500, height = 420*(as.integer(ncol(dd)/3)+1))
par(mfrow=c(as.integer(ncol(dd)/3)+1,3))
for (i in 3: ncol(dd)){
  heading = opts[i]
  plot(dd$POS, dd[,i],
       main=heading,
       ylim = c(min(dd[,-1:-2]),max(dd[,-1:-2])),
       cex.main = 2, cex.axis = 1.5, cex.lab = 1.5,
       col="blue", pch = ".",
       xlab = "Reference position (NT)",ylab = "Read Counts")
  lines(dd$POS, dd[,i], lty = 1, col = "blue", lwd = 1)
  abline(h = 5000, lty = 1, col = "red")
  text(refgenome, x = 800, y = 100, cex = 2)
  }
dev.off()

cov.mean = colMeans(dd[,-1:-2])
cov.mean = as.data.frame(cov.mean)
write.table(cov.mean,
	    file = args[2],
	    quote = F,sep = "\t",row.names = T,col.names = F)
