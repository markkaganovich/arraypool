# R --no-restore --no-save --args pool2inlowcov.lowcov.genocomma Array25M3.Rinput .d34 Array25M4.Rinput Array25M5.Rinput < arrays.R
#setwd('/Users/markkaganovich/arraypool/')
args <- commandArgs(trailingOnly = TRUE)
print(args)

gfile = args[1]
controlsnpfile = args[2]
outputext = args[3]
#expsnpsfiles = args[4:length(args)]
ef = args[4]

g.raw <- read.table(gfile,sep=',', header=F, row.names=1)

csnp.raw <- read.table(controlsnpfile, sep = '\t', header = F, row.names = 1)


#if (length(expsnpsfiles) > 1){
#	for( ef in expsnpsfiles){

		esnp.raw <- read.table(ef, sep = '\t', header = F, row.names = 1)
		snps = intersect(row.names(na.omit(csnp.raw)), row.names(na.omit(esnp.raw)))
		isnps = intersect(snps, row.names(g.raw))
		ind = rownames(g.raw) %in% isnps 
		g = g.raw[ind,]
		g = as.matrix(g)
		ind.c = rownames(csnp.raw) %in% isnps
		csnp = csnp.raw[ind.c,]
		ind.e = rownames(esnp.raw) %in% isnps
		esnp = esnp.raw[ind.e,]

		d <- solve(t(g) %*% g, t(g) %*% (esnp - csnp))

		output = paste(controlsnpfile,ef,outputext, sep='')
		save(d, file=output)

#	}
#}


