snptotal = 7724854
lines = 60

Gceu <-as.matrix(read.table('~/1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcfmatrixoutput',sep="\t", nrows = snptotal))[,-(lines+1)]

lcls.uniform <- rep(1/lines, lines)
S.uniform <- Gceu %*% lcls.uniform
lcls.beta <- solve(t(Gceu) %*% Gceu, t(Gceu) %*% S.uniform)

#generate random test.lcls vectors
simlen <- 100
linesensitivity <- c()
for (j in 1:lines){
    test.lcls <- rep(0,lines)
    test.lcls[j] = 1/100
    for (k in 1:simlen){
	    i = sample(lines,1)
	    if (i != j){
            test.lcls[i] = sample(lines - sum(test.lcls), 1) 
}
}

test.lcls <- test.lcls/lines
S.test <- Gceu %*% test.lcls
test.lcls.dist <- sum((test.lcls - lcls.uniform)^2)

# simulate adding error, for each test.lcls
#test.beta <- solve(t(Gceu) %*% Gceu, t(Gceu) %*% S.test)
#test.residual <- test.beta - test.lcls
#errors <- c()
#for (i in 1:10/50){
std = .05
S.error <- S.test + rnorm(snptotal, mean = 0, sd = std)
test.beta.error <- solve(t(Gceu) %*% Gceu, t(Gceu) %*% S.error)
test.error.residual <- test.beta.error - test.lcls
error <- sum(test.error.residual^2)
linesensitivity[j] <- error
#errors <-c(errors, sum(test.error.residual^2))
#}
}


