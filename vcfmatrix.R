snptotal = 7724854
lines = 60

G <-as.matrix(read.table('~/1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcfmatrixoutput',sep="\t", nrows = snptotal))[,-(lines+1)]

lcls.uniform <- rep(1/lines, lines)
S.uniform <- G %*% lcls.uniform
lcls.beta <- solve(t(G) %*% G, t(G) %*% S.uniform)

#generate random test.lcls vectors
simlen <- 100
simlen2 <- 1000
linesensitivity <- matrix(c(rep(0,180)), nrow=60, ncol=3)
freqs <- c(1/1000,1/100,1/10)
for (j in 1:lines){
    errors <- c()
    for (f in 1:length(freqs)){
        simerrors <- c()
        for (s in 1:simlen){
            test.lcls <- rep(0,lines)
            test.lcls[j] = freqs[f]
            for (k in 1:simlen2){
	            i = sample(lines,1)
	            if (i != j){
                    test.lcls[i] = sample(lines - sum(test.lcls), 1) 
                }
            }

        test.lcls <- test.lcls/lines
        S.test <- G %*% test.lcls
        test.lcls.dist <- sum((test.lcls - lcls.uniform)^2)

# simulate adding error, for each test.lcls
#test.beta <- solve(t(Gceu) %*% Gceu, t(Gceu) %*% S.test)
#test.residual <- test.beta - test.lcls
#errors <- c()
#for (i in 1:10/50){
        std = .05
        S.error <- S.test + rnorm(snptotal, mean = 0, sd = std)
        test.beta.error <- solve(t(G) %*% G, t(G) %*% S.error)
        test.error.residual <- abs(test.beta.error[j] - test.lcls[j])/test.lcls[j]
        simerrors <- c(simerrors, test.error.residual)
    }        
    error <- sum(simerrors) / simlen
    linesensitivity[j , f] = error
    print(j)
    write(linesensitivity, "~/arraypool/test")
#errors <-c(errors, sum(test.error.residual^2))
#}
}   
}

