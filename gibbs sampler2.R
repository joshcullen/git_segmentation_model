rm(list=ls(all=TRUE))
set.seed(1)

setwd('U:\\GIT_models\\git_segmentation_model')
source('gibbs functions.R')
dat=read.csv('fake data.csv',as.is=T)

#priors
alpha=0.01

#useful stuff
max.time=max(dat$time1)
nloc=ncol(dat)-1

#starting values
breakpt=mean(dat$time1)

#number of iterations
ngibbs=10000

#matrix to store results
store.param=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals<- samp.move(breakpt=breakpt,max.time=max.time,dat=dat,
                   alpha=alpha,nloc=nloc)
  breakpt=vals[[1]]
  
  #store results
  store.param[i,]=c(length(vals[[1]]), vals[[2]])
}

length(breakpt)
abline(v=breakpt/nrow(dat),lty=3,col='green')


## MCMC traceplots and diagnostics

#create mcmc object for diags
store.param.mcmc<- as.mcmc(store.param)

#breakpts
traceplot(store.param.mcmc[,1]); abline(h=nseg-1,col='red')

#LML
traceplot(store.param.mcmc[,2])

