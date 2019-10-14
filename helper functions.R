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
grid.summary.table=function(dat,crs){  #dat must already have time.seg assigned
  
  #create grid and extract coords per cell
  grid_5<- raster(extent(min(dat$utmlong), max(dat$utmlong),
                         min(dat$utmlat), max(dat$utmlat)) + 10000)
  res(grid_5)<- 5000
  proj4string(grid_5)<- crs
  grid_5[]<- 0
  grid.cell.locs<- coordinates(grid_5) %>% data.frame()
  names(grid.cell.locs)<- c("grid.x", "grid.y")
  grid.cell.locs$grid.cell<- 1:length(grid_5)
  
  #filter by only cells visited by ID
  dat<- left_join(dat,grid.cell.locs, by="grid.cell")
  
  grid.cell.id<- unique(dat$grid.cell)
  dat.cells=matrix(NA, length(grid.cell.id), 3)
  
  for (i in 1:length(grid.cell.id)) {
    tmp<- dat %>% filter(grid.cell == grid.cell.id[i]) %>% dplyr::select(grid.cell,grid.x,grid.y) %>%
      slice(n=1) %>% as.numeric()
    
    dat.cells[i,]<- tmp
  }
  
  colnames(dat.cells)<- c("grid.cell","grid.x","grid.y")
  
  dat.cells
}
#------------------------------------------------
k.select <- function(){  #to be used by kmeans.cluster()
  k <- readline("What is the value of k?  ")  
  
  k <- as.numeric(unlist(strsplit(k, ",")))
  
  set.seed(1)
  
  dat.kmeans<- kmeans(dat.cells[,-1], k)
  dat.cells<- dat.cells %>% data.frame()
  dat.cells$k.clust<- dat.kmeans$cluster
  
  return(dat.cells)
  
}
#------------------------------------------------
kmeans.cluster=function(dat.cells) {
  store.val=matrix(NA,20,1)
  
  set.seed(1)
  for (i in 1:20) {  #test from 1 to 20 clusters
    tmp=kmeans(dat.cells[,-1], i)
    store.val[i,]=tmp$betweenss / tmp$totss
  }
  
  plot(store.val)
  
  if(interactive()) k.select()
}
#------------------------------------------------
get.summary.stats_kmeans=function(dat){  #dat must have kclust assigned by obs
  kclust=dat$k.clust
  n=length(unique(dat$time.seg))
  res=matrix(0,n,max(dat$k.clust))
  for (i in 1:n){
    ind=dat %>% filter(time.seg==i) %>% group_by(k.clust) %>% count()
    res[i,ind$k.clust]=ind$n #takes count of each cluster within given time segment
  }
  res
}
