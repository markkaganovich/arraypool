# R --no-restore --no-save --args pool2inlowcov.lowcov.geno Array25M3.Rinput Array25M4.Rinput d34 < arrays.R
#setwd('/Users/markkaganovich/arraypool/')
args <- commandArgs(trailingOnly = TRUE)
print(args)

gfile = args[1]
controlsnpfile = args[2]
expsnpsfiles = args[4:length(args)]
outputext = args[3]

g.raw <- read.table(gfile,sep='\t', header=F, row.names=1)
V2 = g.raw$V2
V2 = levels(V2)[V2]
V2 = as.matrix(V2)
row.names(V2) = row.names(g.raw)

csnp.raw <- read.table(controlsnpfile, sep = '\t', header = F, row.names = 1)


if (length(expsnpsfiles) > 1){
	for( ef in expsnpsfiles){
		esnp.raw <- read.table(ef, sep = '\t', header = F, row.names = 1)
		snps = intersect(row.names(na.omit(csnp.raw)), row.names(na.omit(esnp.raw)))
		g = g.m[snps,]
		csnp = csnp.raw[snps,]
		esnp = esnp.raw[snps,]

		d <- solve(t(g) %*% g, t(g) %*% (esnp - csnp))

		output = paste(controlsnpfile,ef,outputext, sep='')
		save(d, output)

	}
}




#print(gfile)
#print(controlsnpfile)
#print(expsnpsfiles[1])
#print(outputext)