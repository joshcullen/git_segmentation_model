assign.time.seg=function(dat){
  tmp=which(unique(dat$id) == identity)
  breakpt<- dat.res$brkpts[tmp,-1] %>% discard(is.na) %>% as.numeric(.[1,])
  breakpt1=c(0,breakpt,Inf)
  n=length(breakpt1)
  res=matrix(0,nrow(dat),1)
  for (i in 2:n){
    ind=which(breakpt1[i-1]<dat$time1 & dat$time1<breakpt1[i])
    res[ind,]=i-1
  }
  dat$time.seg<- as.vector(res)
  dat
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
#------------------------------------------------
traceplot=function(data, type) {  #create traceplots for nbrks or LML for all IDs
  for (i in 1:length(identity)) {
    par(ask=TRUE)
    plot(x=1:ngibbs, y=data[i,-1], type = "l", xlab = "Iteration",
         ylab = ifelse(type == "nbrks", "# of Breakpoints", "Log Marginal Likelihood"),
         main = paste("ID",data[i,1]))
  }
  on.exit(par(ask = FALSE))
}
#------------------------------------------------
plot.heatmap=function(data) {
  
  #re-define loc.id based only on those visited by this individual
  uni.loc=unique(data$loc.id)
  aux=data.frame(loc.id=uni.loc,loc.id1=1:length(uni.loc))
  dat1=merge(data,aux,all=T)
  dat1$loc.id=dat1$loc.id1
  dat=dat1[order(dat1$time1),c('loc.id','time1')]
  
  nloc<- length(uni.loc)
  nobs<- nrow(data)
  obs<- matrix(0, nobs, nloc)
  
  for (i in 1:nrow(data)) {
    obs[i, dat$loc.id[i]]<- 1
  }
  
  obs<- data.frame(obs)
  names(obs)<- 1:nloc
  obs.long<- obs %>% gather(key, value) %>% mutate(time=rep(1:nobs, times=nloc))
  obs.long$key<- as.numeric(obs.long$key)
  
  tmp=which(unique(data$id) == identity)
  breakpt<- dat.res$brkpts[tmp,-1] %>% discard(is.na) %>% t() %>% data.frame()
  names(breakpt)<- "breaks"
  
  
  print(
    ggplot(obs.long, aes(x=time, y=key, fill=value)) +
      geom_tile() +
      scale_fill_viridis_c(guide=F) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      geom_vline(data = breakpt, aes(xintercept = breaks), color = viridis(n=9)[7], size = 0.35) +
      labs(x = "Observations", y = "Grid Cell") +
      theme_bw() +
      theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16)) +
      labs(title = paste("ID", unique(data$id)))
  )
  
  
}
#------------------------------------------------
heatmap=function(data) {
  par(ask = TRUE)
  map(data, plot.heatmap)
  par(ask = FALSE)
}
