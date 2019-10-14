ind.func=function(x) {which(x==1)}
#------------------------------------------------


get.summary.stats=function(breakpt,dat,nloc){
  col.time1=which(colnames(dat)=='time1')
  breakpt1=c(0,breakpt,Inf)
  n=length(breakpt1)
  res=matrix(NA,n-1,nloc)
  for (i in 2:n){
    ind=which(breakpt1[i-1]<dat$time1 & dat$time1<breakpt1[i])
    tmp=dat[ind,-col.time1]
    
    ind.loc<- matrix(NA, nrow(tmp), 1)
    ind.loc<- apply(tmp, 1, FUN = ind.func)  #creates vector of locations occupied per observation
    
    tab<- tabulate(ind.loc)  #ensures that there are values for all nloc
    if (length(tab) < nloc) {
      tab<- c(tab,rep(0,nloc-length(tab)))
    }
    
    res[i-1,]=tab #takes count of the each location within given time segment
  }
  res
}
#---------------------------------------------
log.marg.likel=function(alpha,summary.stats,nloc){
  #get ratio
  lnum=rowSums(lgamma(alpha+summary.stats))
  lden=lgamma(nloc*alpha+rowSums(summary.stats))
  p2=sum(lnum)-sum(lden)
  p1=nrow(summary.stats)*(lgamma(nloc*alpha)-nloc*lgamma(alpha))
  p1+p2
}
#---------------------------------------------
samp.move=function(breakpt,max.time,dat,alpha,nloc){
  breakpt.old=breakpt
  p=length(breakpt)
  rand1=runif(1)	
  p0=1
  new.brk=runif(1,min=0,max=max.time)
  
  if (p == 1) {
    #birth
    if (rand1 < 1/2){
      breakpt.new=sort(c(breakpt.old,new.brk))
      p0=2/3 #death prob 2 -> 1 is (1/3) and birth prob 1 -> 2 is 1/2. 
    }
    #swap
    if (rand1 > 1/2) breakpt.new=new.brk
  }
  if (p > 1) {
    #birth
    if (rand1 < 1/3) {
      breakpt.new=sort(c(breakpt.old,new.brk))
    }
    #death
    if (rand1 > 1/3 & rand1 < 2/3) {
      ind=sample(1:length(breakpt.old),size=1)
      breakpt.new=breakpt.old[-ind]
      if (p==2) p0=3/2 #birth prob from 1 -> 2 is 1/2 and death prob from 2 -> 1 is 1/3
    }
    #swap
    if (rand1 > 2/3) {
      ind=sample(1:length(breakpt.old),size=1)
      breakpt.new=sort(c(breakpt.old[-ind],new.brk))
    }
  }
  
  #get sufficient statistics
  stats.old=get.summary.stats(breakpt=breakpt.old,dat=dat,nloc=nloc)
  stats.new=get.summary.stats(breakpt=breakpt.new,dat=dat,nloc=nloc)
  
  #get marginal loglikel
  pold=log.marg.likel(alpha=alpha,summary.stats=stats.old,nloc=nloc)
  pnew=log.marg.likel(alpha=alpha,summary.stats=stats.new,nloc=nloc)+log(p0)
  prob=exp(pnew-pold)
  rand2=runif(1)
  
  if (rand2<prob) return(list(breakpt.new, pnew))
  return(list(breakpt.old, pold))
  
}

