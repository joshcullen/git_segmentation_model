assign.time.seg<- function(obs, breakpts, dat) {
  seg.n=apply(obs,1,sum)
  time.seg=rep(1:(length(breakpts)+1), seg.n)
  dat$time.seg<- time.seg
  
  dat
}
#------------------------------------------------
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
#------------------------------------------------
create.grid=function(dat,crs,extent,res,buffer) {
  grid<- raster(extent + buffer)
  res(grid)<- res
  proj4string(grid)<- crs
  grid[]<- 0
  
  grid
}
#------------------------------------------------
grid.summary.table=function(dat,crs,extent,res,buffer){  #dat must already have time.seg assigned
  
  #create grid and extract coords per cell
  grid<- create.grid(dat=dat, crs=crs, extent=extent, res=res, buffer=buffer)
  grid.cell.locs<- coordinates(grid) %>% data.frame()
  names(grid.cell.locs)<- c("x", "y")
  grid.cell.locs$grid.cell<- 1:length(grid)
  grid.coord<- grid.cell.locs[grid.cell.locs$grid.cell %in% dat$grid.cell,]
  
  
  grid.coord
}
#------------------------------------------------
df.to.list=function(dat) {  #only for id as col in dat
    id<- unique(dat$id)
    n=length(id)
    dat.list<- vector("list", n)
    names(dat.list)<- id
    
    for (i in 1:length(id)) {
      dat.list[[i]]<- dat[dat$id==id[i],]
    }
    dat.list
}