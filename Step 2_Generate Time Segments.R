#Analyze Snail Kite Data with Bayesian Partitioning Model

set.seed(1)

library(dplyr)
library(lubridate)
library(ggplot2)
library(coda)
library(raster)

setwd("~/Documents/Snail Kite Project/Cat-Dirichlet Gibbs Sampler")
source('gibbs functions2.R')

setwd("~/Documents/Snail Kite Project/Data")
dat<- read.csv("Snail Kite Gridded Data.csv", header = T, sep = ",")

#Update class and/or values of vars
dat$ESTtime<- as_datetime(dat$ESTtime)

############################################
#### Create Grid Occupancy Matrix by ID ####
############################################

#Create grid for indexing
grid_5<- raster(extent(min(dat$utmlong), max(dat$utmlong),
                         min(dat$utmlat), max(dat$utmlat)) + 10000)
res(grid_5)<- 5000
proj4string(grid_5)<- CRS("+init=epsg:32617")
grid_5[]<- 0

#Create occupancy matrices
obs<- matrix(0, nrow(dat), length(grid_5))
for (i in 1:nrow(dat)){
  obs[i,dat$grid.cell[i]]=1
}

colnames(obs)=paste0('grid.cell_',1:length(grid_5))
obs<- obs[,apply(obs,2,sum) != 0]
obs<- cbind(dat$id, obs)
colnames(obs)[1]<- "id"

write.csv(obs, "Occupancy Matrix for all Obs and Locs.csv", row.names = F)


obs1=data.frame(obs) %>% filter(id == 1) %>% dplyr::select(-id)
obs12=data.frame(obs) %>% filter(id == 12) %>% dplyr::select(-id)
obs19=data.frame(obs) %>% filter(id == 19) %>% dplyr::select(-id)
obs27=data.frame(obs) %>% filter(id == 27) %>% dplyr::select(-id)

image(as.matrix(obs1)); title(x="Observations (n=9134)", y="Locations (n=378)")
image(as.matrix(obs12)); title(x="Observations (n=4922)", y="Locations (n=378)")
image(as.matrix(obs19)); title(x="Observations (n=3068)", y="Locations (n=378)")
image(as.matrix(obs27)); title(x="Observations (n=1182)", y="Locations (n=378)")

obs1$time1=1:nrow(obs1)
obs12$time1=1:nrow(obs12)
obs19$time1=1:nrow(obs19)
obs27$time1=1:nrow(obs27)




#################################
#### Run Gibbs Sampler by ID ####
#################################

nloc=ncol(obs)-1


### ID 1

#priors
alpha=0.01

#useful stuff
max.time=max(obs1$time1)

#starting values
breakpt=mean(obs1$time1)

#number of iterations
ngibbs=10000

#matrix to store results
store.param1=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals<- samp.move(breakpt=breakpt,max.time=max.time,dat=obs1,
                   alpha=alpha,nloc=nloc)
  breakpt=vals[[1]]
  
  #store results
  store.param1[i,]=c(length(vals[[1]]), vals[[2]])
}
###Takes ~ 1.64 hrs to run for 10000 iterations; identified 48 breakpoints
##LML of final model is -5987.567


length(breakpt)
write.csv(breakpt, "ID1 Breakpoints (5 km).csv", row.names = F)

setwd("~/Documents/Snail Kite Project/Data/Figures")
jpeg("Time Segments and Breakpoints for ID1.jpeg", width = 6, height = 4, units = "in", res = 400)
image(as.matrix(obs1[,-379])); title(x="Observations (n=9134)", y="Locations (n=378)")
abline(v=breakpt/nrow(obs1),lty=3,col='blue')
dev.off()


## MCMC traceplots and diagnostics

#create mcmc object for diags
store.param1.mcmc<- as.mcmc(store.param1)

#breakpts
jpeg("Traceplot of Breakpoints for ID1.jpeg", width = 6, height = 4, units = "in", res = 400)
traceplot(store.param1.mcmc[,1])
title(y="# of Breakpoints")
dev.off()

#LML
jpeg("Traceplot of LML for ID1.jpeg", width = 6, height = 4, units = "in", res = 400)
traceplot(store.param1.mcmc[,2])
title(y="Log Marginal Likelihood")
dev.off()











### ID 12

#priors
alpha=0.01

#useful stuff
max.time=max(obs12$time1)

#starting values
breakpt=mean(obs12$time1)

#number of iterations
ngibbs=10000

#matrix to store results
store.param12=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals<- samp.move(breakpt=breakpt,max.time=max.time,dat=obs12,
                   alpha=alpha,nloc=nloc)
  breakpt=vals[[1]]
  
  #store results
  store.param12[i,]=c(length(vals[[1]]), vals[[2]])
}
###Takes ~ 1.13 hrs to run for 10000 iterations; identified 47 breakpoints
##LML of final model is -3484.219


length(breakpt)
write.csv(breakpt, "ID12 Breakpoints (5 km).csv", row.names = F)

jpeg("Time Segments and Breakpoints for ID12.jpeg", width = 6, height = 4, units = "in", res = 400)
image(as.matrix(obs12[,-379])); title(x="Observations (n=4922)", y="Locations (n=378)")
abline(v=breakpt/nrow(obs12),lty=3,col='blue')
dev.off()


## MCMC traceplots and diagnostics

#create mcmc object for diags
store.param12.mcmc<- as.mcmc(store.param12)

#breakpts
jpeg("Traceplot of Breakpoints for ID12.jpeg", width = 6, height = 4, units = "in", res = 400)
traceplot(store.param12.mcmc[,1])
title(y="# of Breakpoints")
dev.off()

#LML
jpeg("Traceplot of LML for ID12.jpeg", width = 6, height = 4, units = "in", res = 400)
traceplot(store.param12.mcmc[,2])
title(y="Log Marginal Likelihood")
dev.off()











### ID 19

#priors
alpha=0.01

#useful stuff
max.time=max(obs19$time1)

#starting values
breakpt=mean(obs19$time1)

#number of iterations
ngibbs=10000

#matrix to store results
store.param19=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals<- samp.move(breakpt=breakpt,max.time=max.time,dat=obs19,
                   alpha=alpha,nloc=nloc)
  breakpt=vals[[1]]
  
  #store results
  store.param19[i,]=c(length(vals[[1]]), vals[[2]])
}
###Takes ~ 0.89 hrs to run for 10000 iterations; identified 35 breakpoints
##LML of final model is -2604.251


length(breakpt)
write.csv(breakpt, "ID19 Breakpoints (5 km).csv", row.names = F)

jpeg("Time Segments and Breakpoints for ID19.jpeg", width = 6, height = 4, units = "in", res = 400)
image(as.matrix(obs19[,-379])); title(x="Observations (n=3068)", y="Locations (n=378)")
abline(v=breakpt/nrow(obs19),lty=3,col='blue')
dev.off()


## MCMC traceplots and diagnostics

#create mcmc object for diags
store.param19.mcmc<- as.mcmc(store.param19)

#breakpts
jpeg("Traceplot of Breakpoints for ID19.jpeg", width = 6, height = 4, units = "in", res = 400)
traceplot(store.param19.mcmc[,1])
title(y="# of Breakpoints")
dev.off()

#LML
jpeg("Traceplot of LML for ID19.jpeg", width = 6, height = 4, units = "in", res = 400)
traceplot(store.param19.mcmc[,2])
title(y="Log Marginal Likelihood")
dev.off()











### ID 27

#priors
alpha=0.01

#useful stuff
max.time=max(obs27$time1)

#starting values
breakpt=mean(obs27$time1)

#number of iterations
ngibbs=10000

#matrix to store results
store.param27=matrix(NA,ngibbs,2)

for (i in 1:ngibbs){
  print(i)
  vals<- samp.move(breakpt=breakpt,max.time=max.time,dat=obs27,
                   alpha=alpha,nloc=nloc)
  breakpt=vals[[1]]
  
  #store results
  store.param27[i,]=c(length(vals[[1]]), vals[[2]])
}
###Takes ~ 0.28 hrs to run for 10000 iterations; identified 7 breakpoints
##LML of final model is -984.5128


length(breakpt)
write.csv(breakpt, "ID27 Breakpoints (5 km).csv", row.names = F)

jpeg("Time Segments and Breakpoints for ID27.jpeg", width = 6, height = 4, units = "in", res = 400)
image(as.matrix(obs27[,-379])); title(x="Observations (n=1182)", y="Locations (n=378)")
abline(v=breakpt/nrow(obs27),lty=3,col='blue')
dev.off()


## MCMC traceplots and diagnostics

#create mcmc object for diags
store.param27.mcmc<- as.mcmc(store.param27)

#breakpts
jpeg("Traceplot of Breakpoints for ID27.jpeg", width = 6, height = 4, units = "in", res = 400)
traceplot(store.param27.mcmc[,1])
title(y="# of Breakpoints")
dev.off()

#LML
jpeg("Traceplot of LML for ID27.jpeg", width = 6, height = 4, units = "in", res = 400)
traceplot(store.param27.mcmc[,2])
title(y="Log Marginal Likelihood")
dev.off()



#####################################################
#### Import Original Model Input and Breakpoints ####
#####################################################

setwd("~/Documents/Snail Kite Project/Data")
obs<- read.csv("Occupancy Matrix for all Obs and Locs.csv", header = T, sep = ",")
obs1.breakpts<- read.csv("ID1 Breakpoints (5 km).csv", header = T, sep = ",")
obs1.breakpts=obs1.breakpts[,1]
obs12.breakpts<- read.csv("ID12 Breakpoints (5 km).csv", header = T, sep = ",")
obs12.breakpts=obs12.breakpts[,1]
obs19.breakpts<- read.csv("ID19 Breakpoints (5 km).csv", header = T, sep = ",")
obs19.breakpts=obs19.breakpts[,1]
obs27.breakpts<- read.csv("ID27 Breakpoints (5 km).csv", header = T, sep = ",")
obs27.breakpts=obs27.breakpts[,1]


obs1=data.frame(obs) %>% filter(id == 1) %>% dplyr::select(-id)
obs12=data.frame(obs) %>% filter(id == 12) %>% dplyr::select(-id)
obs19=data.frame(obs) %>% filter(id == 19) %>% dplyr::select(-id)
obs27=data.frame(obs) %>% filter(id == 27) %>% dplyr::select(-id)

obs1$time1=1:nrow(obs1)
obs12$time1=1:nrow(obs12)
obs19$time1=1:nrow(obs19)
obs27$time1=1:nrow(obs27)

############################################
#### Summarize Results From First Model ####
############################################

nloc=ncol(obs)-1 #remove time1


obs1.seg=get.summary.stats(obs1.breakpts,obs1,nloc)
write.csv(obs1.seg, "ID1 Seg x Loc.csv", row.names = F)
obs12.seg=get.summary.stats(obs12.breakpts,obs12,nloc)
write.csv(obs12.seg, "ID12 Seg x Loc.csv", row.names = F)
obs19.seg=get.summary.stats(obs19.breakpts,obs19,nloc)
write.csv(obs19.seg, "ID19 Seg x Loc.csv", row.names = F)
obs27.seg=get.summary.stats(obs27.breakpts,obs27,nloc)
write.csv(obs27.seg, "ID27 Seg x Loc.csv", row.names = F)


################################################
#### Summarize Results With Site Clustering ####
################################################

dat1=dat %>% filter(id==1)
dat1<- assign.time.seg(obs1.seg, obs1.breakpts, dat1)
dat12=dat %>% filter(id==12)
dat12<- assign.time.seg(obs12.seg, obs12.breakpts, dat12)
dat19=dat %>% filter(id==19)
dat19<- assign.time.seg(obs19.seg, obs19.breakpts, dat19)
dat27=dat %>% filter(id==27)
dat27<- assign.time.seg(obs27.seg, obs27.breakpts, dat27)

dat.list<- list(`1`=dat1, `12`=dat12, `19`=dat19, `27`=dat27)

#for use by clustering of segments with TSB prior
#dat must already have time.seg column

kmeans.list<- vector("list", length(dat.list))

for (i in 1:length(dat.list)) {
dat.cells<- grid.summary.table(dat = dat.list[[i]], crs = CRS('+init=epsg:32617'))  #not working for IDs 12-27; check this function
dat.cells<- kmeans.cluster(dat.cells)
dat<- left_join(dat.list[[i]], dat.cells[,-c(2:3)], by = "grid.cell")
obs.kmeans<- get.summary.stats_kmeans(dat)

kmeans.list[[i]]<- obs.kmeans
}
