#setwd('/Users/markkaganovich/arraypool/')
args <- commandArgs(trailingOnly = TRUE)
print(args)

gfile = args[1]
controlsnpfile = args[2]
expsnpsfile = args[3]

#gfile = '30/poolgenotype.Rinput'
#controlsnpfile = '30/newtheta1'
#expsnpsfile = '30/newtheta3'

g.raw <- read.table(gfile,sep='\t', header=F, row.names=1)
csnp.raw <- read.table(controlsnpfile, sep = '\t', header = F, row.names = 1)
esnp.raw <- read.table(expsnpsfile, sep = '\t', header = F, row.names = 1)

snps = intersect(row.names(na.omit(csnp.raw)), row.names(na.omit(esnp.cor)))

g = as.matrix(g.raw[snps,])
csnp = csnp.raw[snps,]
esnp = esnp.raw[snps,]

d <- solve(t(g) %*% g, t(g) %*% (esnp - csnp))

