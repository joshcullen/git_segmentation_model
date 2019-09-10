rm(list=ls(all=TRUE))
library('MCMCpack')
set.seed(1)

nobs=1000
nseg=10
tmp=runif(nseg)
prob=tmp/sum(tmp); prob
partition=rmultinom(1,size=nobs,prob=prob)
seg.index=rep(1:nseg,times=partition)

nloc=100
prob=rdirichlet(nseg,rep(0.01,nloc))

obs=matrix(NA,nobs,nloc)
for (i in 1:nobs){
  obs[i,]=rmultinom(1,size=1,prob=prob[seg.index[i],])
}
image(obs)
abline(v=cumsum(partition)/nobs)

colnames(obs)=paste0('loc',1:nloc)
obs1=data.frame(obs)
obs1$time1=1:nobs

setwd('U:\\GIT_models\\git_segmentation_model')
write.csv(obs1,'fake data.csv',row.names=F)
