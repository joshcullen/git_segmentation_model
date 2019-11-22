rm(list=ls(all=TRUE))
set.seed(1)

#get functions
setwd('U:\\GIT_models\\git_segmentation_model')
source('gibbs functions.R')

#get data
setwd('U:\\poli\\josh\\tsegmen_loc')
k=read.csv('Snail Kite Data_condensed.csv')
uni.id=unique(k$id)
n.id=length(uni.id)

#priors
alpha=0.01

#to store results
res.gibbs=matrix(NA,n.id,100)
res.time=matrix(NA,n.id,3)

ngibbs=10000
oo=1

for (id in uni.id){
  print(id)
  start_time <- Sys.time() #start time
  
  #subset data for individual id
  cond=k$id==id
  dat=k[cond,c('loc.id','time1')]
  
  #re-define loc.id based only on those visited by this individual
  uni.loc=unique(dat$loc.id)
  aux=data.frame(loc.id=uni.loc,loc.id1=1:length(uni.loc))
  dat1=merge(dat,aux,all=T); dim(dat); dim(dat1)
  dat1$loc.id=dat1$loc.id1
  dat=dat1[order(dat1$time1),c('loc.id','time1')]
  
  #useful stuff
  max.time=max(dat$time1)
  nloc=max(dat$loc.id)
  
  #starting values
  breakpt=mean(dat$time1)

  for (i in 1:ngibbs){
    print(i)
    breakpt=samp.move(breakpt=breakpt,max.time=max.time,dat=dat,
                      alpha=alpha,nloc=nloc)   
  }
  end_time <- Sys.time()
  
  tmp=c(id,breakpt)
  res.gibbs[oo,1:length(tmp)]=tmp
  tmp=end_time - start_time
  res.time[oo,]=c(id,tmp,units(tmp))
  oo=oo+1
} 
colnames(res.time)=c('id','elapsed time','units')
res.time[,2]=round(as.numeric(res.time[,2]),3)
res.time




