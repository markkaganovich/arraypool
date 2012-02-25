snptotal = 7724854

Gceu <- as.matrix(read.table('~/1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcfmatrixoutput', sep="\t", nrows = snptotal))[,-61]
Sceu <- Gceu/60

lcls <- rep(1/60, 60)

model <- glm(lcls ~ Sceu)

# check the math, make a simulated genotype matrix that would be measure in a pool, check that the linear regression gives us back the inputed lcl frequency matrix

testlcls <- lcls
testlcls[1] = 1/120
testlcls[2] = 1/120
testlcls[3] = 1/30

testS <- Gceu %*% testlcls

model <- lm(testS ~ Gceu)

beta.hat = solve(t(Gceu) %*% Gceu, t(Gceu) %*% testS)
beta.hat.new = solve(t(Gceu) %*% Gceu, t(Gceu) %*% S.error)


#generate random L vectors
simlen <- 100
len <- 60
L <-rep(0, len)
for (k in 1:simlen){
	i = sample(len,1)
	L[i] = sample(len - sum(L), 1) 
}
L <- L/len

S <- Gceu %*% L


# simulate adding error
S.error <- S
for (i in 1:60){
    error<-rnorm(snptotal, mean = 0, sd=.05)
    S.error[,i] <- S[,i] + error
}



errorsimlen = 10000
S.error <- S
for (k in 1:errorsimlen){
e <- sample(100,1)/100
i = sample(snptotal,1)
#S.error[i] = e* S[i]
S.error[i] = 10 * S[i]
}
 


