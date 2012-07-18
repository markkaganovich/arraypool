#G <- as.matrix(read.table('./Gnotrios', sep=',', nrows = 473822))

topten.score <- rep(0:10)
topten.id <- c()
topten.lcls <- c()

sim=0
while(sim < 100){
lcls<-c()
while(sum(lcls) < 90 | sum(lcls) > 110){
  lcls.id <- sample(1:179,sample(2:179,1))
  lcls <- sample(1:100, length(lcls.id), replace=TRUE)
}

S.sim <- (G[, lcls.id] %*% (lcls/100))/(length(lcls) *2)
v <- as.matrix(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
S.goodsnps <- apply(v, 1, function(x) length(which(S.sim > (x-.01) & S.sim < (x+.01))))
index = S.goodsnps > topten.score

#keep lcls that score the highest in each dimension
topten.score[index] = S.goodsnps[index]
topten.id[index] = as.data.frame(lcls.id)
topten.lcls[index] = as.data.frame(lcls)

sim = sim +1
}

sink('./simulationoutput')