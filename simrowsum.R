gfile = 'poolgenotype3.Rinput'
ef = '25M1.3.Rinput'
controlsnpfile = '25M1.1.Rinput'

spiked <- c(8,4,8,.25,32,8,8,8,8,4,4,4,4,32,16,2,2,2,2,2,32,16,32,1,16,4,1,1)
spiked1 = spiked/sum(spiked)

g.raw <- read.table(gfile,sep=',', header=F, row.names=1)
csnp.raw <- read.table(controlsnpfile, sep = '\t', header = F, row.names = 1)
esnp.raw <- read.table(ef, sep = '\t', header = F, row.names = 1)
snps = intersect(row.names(na.omit(csnp.raw)), row.names(na.omit(esnp.raw)))
isnps = intersect(snps, row.names(g.raw))

scores = c()

for (size in seq(1,59)){
	s = names(which(rowSums(g.raw) >= size))
	isnps.s = intersect(isnps, s)
	ind = rownames(g.raw) %in% isnps.s
	ind.c = rownames(csnp.raw) %in% isnps.s
	ind.e = rownames(esnp.raw) %in% isnps.s
	g = g.raw[ind,]
	g = as.matrix(g)
	csnp = csnp.raw[ind.c,]
	esnp = esnp.raw[ind.e,]
	d <- solve(t(g) %*% g, t(g) %*% (esnp - csnp))

	reald = (d+1/28)/sum(d+1/28)
	scores[size] = sum((reald - spiked1)^2)
	print(size)
	print(sum((reald - spiked1)^2))


}