gibbs.time.seg=function(dat,ngibbs) {
set.seed(1)

#priors
alpha=0.01

#useful stuff
max.time=max(dat$time1)
nloc=ncol(dat)-1

#starting values
breakpt=mean(dat$time1)

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
list(breakpt=breakpt, store.param=store.param)
}
