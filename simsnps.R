


spiked <- c(8,4,8,.25,32,8,8,8,8,4,4,4,4,32,16,2,2,2,2,2,32,16,32,1,16,4,1,1)
spiked1 = spiked/sum(spiked)

g.raw <- read.table(gfile,sep=',', header=F, row.names=1)
csnp.raw <- read.table(controlsnpfile, sep = '\t', header = F, row.names = 1)
esnp.raw <- read.table(ef, sep = '\t', header = F, row.names = 1)
snps = intersect(row.names(na.omit(csnp.raw)), row.names(na.omit(esnp.raw)))
isnps = intersect(snps, row.names(g.raw))


meanscores = c()

for (size in range(1000:1000:10000)){
	
	scores = c()
	for (sim in range(1:100)){

	isnps.s = isnps[sample(1:length(isnps), size, ,replacement=FALSE)]

	ind = rownames(g.raw) %in% isnps 
	g = g.raw[ind.s,]
	g = as.matrix(g)		
	ind.c = rownames(csnp.raw) %in% isnps
	csnp = csnp.raw[ind.c,]
	ind.e = rownames(esnp.raw) %in% isnps
	esnp = esnp.raw[ind.e,]

	d <- solve(t(g) %*% g, t(g) %*% (esnp - csnp))

	reald = (d+1/28)/sum(d)
	scores[sim] = reald^2 - d^2
}

print(mean(scores))

meanscores[length(meanscores)+1] = mean(scores)
}
