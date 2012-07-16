G <- as.matrix(read.table('./Gnotrios', sep=',', nrows = 473822))

topten.score <- rep(0:9)
topten.inds <- c()
topten.id <- c()
topten.lcls <- c()

sim=0
while(sim < 10000){
lcls<-c()
while(sum(lcls) < 90 | sum(lcls) > 110){
  lcls.id <- sample(1:179,sample(2:179))
  lcls <- sample(1:100, length(lcls.id), replace=TRUE)
}

S.sim <- (G[, lcls.id] %*% (lcls/100))/(length(lcls) *2)

S.goodsnps <- which((S.sim <.11 & S.sim >.09) | (S.sim < .21 &  S.sim > .19) | (S.sim < .31 & S.sim > .29) |
  (S.sim < .41 & S.sim > .39) | (S.sim < .51 & S.sim > .49) | (S.sim < .61 & S.sim > .59) |
  (S.sim < .71 & S.sim > .69) | (S.sim < .81 & S.sim > .79) | (S.sim < .91 & S.sim > .89)
)

#keep top 10 of these numbers
# also calc dist to maximize dist among the top 10 or don't include those with indeces already covered
if (length(S.goodsnps) > min(topten.score)){
  index = which(topten.score == min(topten.score))
  topten.score[index] = length(S.goodsnps)
  topten.inds[index] = as.data.frame(S.goodsnps)
  topten.id[index] = as.data.frame(lcls.id)
  topten.lcls[index] = as.data.frame(lcls)
}
 
sim = sim +1
}
  