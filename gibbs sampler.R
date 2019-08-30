rm(list=ls(all=TRUE))
set.seed(1)

setwd('U:\\GIT_models\\git_segmentation_model')
source('gibbs functions.R')
dat=read.csv('fake data.csv',as.is=T)

#priors
tau2=100
mu0=0

#useful stuff
max.time=max(dat$time1)

#starting values
breakpt=mean(dat$time1)

ngibbs=10000
for (i in 1:ngibbs){
  print(i)
  breakpt=samp.move(breakpt=breakpt,max.time=max.time,dat=dat,tau2=tau2,mu0=mu0)   
}
length(breakpt)
abline(v=breakpt,col='grey')