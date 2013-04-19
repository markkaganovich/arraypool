
gfile = 'poolgenotype3.Rinput'
ef = '25M1.3.Rinput'
controlsnpfile = '25M1.1.Rinput'

spiked <- c(8,4,8,.25,328,8,8,8,4,4,4,4,32,16,2,2,2,2,2,32,16,32,1,16,4,1,1)
spiked1 = spiked/sum(spiked)

g.raw <- read.table(gfile,sep=',', header=F, row.names=1)
csnp.raw <- read.table(controlsnpfile, sep = '\t', header = F, row.names = 1)
esnp.raw <- read.table(ef, sep = '\t', header = F, row.names = 1)
snps = intersect(row.names(na.omit(csnp.raw)), row.names(na.omit(esnp.raw)))
isnps = intersect(snps, row.names(g.raw))


meanscores = c()

for (size in seq(1000,10000, by = 1000)){
	
	scores = c()
	for (sim in seq(1,100)){

	isnps.s = isnps[sample(length(isnps), size, ,replace=FALSE)]

	ind = rownames(g.raw) %in% isnps.s 
	g = g.raw[ind,]
	g = as.matrix(g)		
	ind.c = rownames(csnp.raw) %in% isnps.s
	csnp = csnp.raw[ind.c,]
	ind.e = rownames(esnp.raw) %in% isnps.s
	esnp = esnp.raw[ind.e,]

	d <- solve(t(g) %*% g, t(g) %*% (esnp - csnp))

	reald = (d+1/28)/sum(d+1/28)
	scores[sim] = sum((reald - spiked1)^2)
}

print(mean(scores))

meanscores[length(meanscores)+1] = mean(scores)
}

save(meanscores, file='simoutputmeanscores1000')

for (size in seq(10000,100000, by = 10000)){
	
	scores = c()
	for (sim in seq(1,100)){

	isnps.s = isnps[sample(length(isnps), size, ,replace=FALSE)]

	ind = rownames(g.raw) %in% isnps.s 
	g = g.raw[ind,]
	g = as.matrix(g)		
	ind.c = rownames(csnp.raw) %in% isnps.s
	csnp = csnp.raw[ind.c,]
	ind.e = rownames(esnp.raw) %in% isnps.s
	esnp = esnp.raw[ind.e,]

	d <- solve(t(g) %*% g, t(g) %*% (esnp - csnp))

	reald = (d+1/28)/sum(d+1/28)
	scores[sim] = sum((reald - spiked1)^2)
}

print(mean(scores))

meanscores[length(meanscores)+1] = mean(scores)
}

save(meanscores, file='simoutputmeanscores10000')


