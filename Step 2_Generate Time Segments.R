#Analyze Snail Kite Data with Bayesian Partitioning Model

library(dplyr)
library(ggplot2)
library(raster)

source('gibbs functions2.R')
source('helper functions.R')
source('gibbs sampler2.R')

dat<- read.csv("Snail Kite Gridded Data_large.csv", header = T, sep = ",")
dat.list<- df.to.list(dat = dat)



############################################
#### Create Grid Occupancy Matrix by ID ####
############################################

#Create grid for indexing
utm.crs<- CRS("+init=epsg:32617")
extent<- extent(min(dat$utmlong), max(dat$utmlong), min(dat$utmlat), max(dat$utmlat))
res<- 5000
buffer<- 10000
grid_5<- create.grid(dat=dat, crs=utm.crs, extent=extent, res=res, buffer=buffer)
#Create occupancy matrices
obs<- matrix(0, nrow(dat), length(grid_5))
for (i in 1:nrow(dat)){
  obs[i,dat$grid.cell[i]]=1
}

colnames(obs)=paste0('grid.cell_',1:length(grid_5))
obs<- obs[,apply(obs,2,sum) != 0]
obs<- cbind(dat$id, obs)
colnames(obs)[1]<- "id"

#write.csv(obs, "Occupancy Matrix for all Obs and Locs.csv", row.names = F)


obs1=data.frame(obs) %>% filter(id == 1) %>% dplyr::select(-id) %>% mutate(time1=1:nrow(.))
obs12=data.frame(obs) %>% filter(id == 12) %>% dplyr::select(-id) %>% mutate(time1=1:nrow(.))
obs19=data.frame(obs) %>% filter(id == 19) %>% dplyr::select(-id) %>% mutate(time1=1:nrow(.))
obs27=data.frame(obs) %>% filter(id == 27) %>% dplyr::select(-id) %>% mutate(time1=1:nrow(.))

nloc=ncol(obs)-1




### Alternative Format for Data ###

dat.list<- lapply(dat.list, function(x) x %>% mutate(time1 = 1:nrow(x)))  #add row for obs number

dat.ind<- map_dfr(dat.list, `[`) %>% dplyr::select(id, grid.cell, time1)  #create DF
names(dat.ind)[2]<- "loc.id"
dat.ind$loc.id<- dat.ind$loc.id %>% factor()
levels(dat.ind$loc.id)<- 1:length(unique(dat.ind$loc.id))  #change from raw to modified cell ID

# write.csv(dat.ind, "Snail Kite Data_condensed.csv", row.names = F)





#################################
#### Run Gibbs Sampler by ID ####
#################################

### ID 1

dat1.res<- gibbs.time.seg(dat.list$`1`, 10000)
###Takes ~ 1.64 hrs to run for 10000 iterations; identified 48 breakpoints
##LML of final model is -5987.567


length(dat1.res$breakpt)
#write.csv(dat1.res$breakpt, "ID1 Breakpoints (5 km).csv", row.names = F)

image(as.matrix(obs1[,-379])); title(x="Observations (n=9134)", y="Locations (n=378)")
abline(v=breakpt/nrow(obs1),lty=3,col='blue')


## traceplots

#breakpts
plot(dat1.res$store.param[[1]], type='l')
title(y="# of Breakpoints")

#LML
plot(dat1.res$store.param[[2]], type='l')
title(y="Log Marginal Likelihood")





### ID 12

dat12.res<- gibbs.time.seg(dat.list$`12`, 10000)
###Takes ~ 1.13 hrs to run for 10000 iterations; identified 47 breakpoints
##LML of final model is -3484.219


length(dat12.res$breakpt)
#write.csv(dat12.res$breakpt, "ID12 Breakpoints (5 km).csv", row.names = F)

image(as.matrix(obs12[,-379])); title(x="Observations (n=4922)", y="Locations (n=378)")
abline(v=breakpt/nrow(obs12),lty=3,col='blue')


## traceplots

#breakpts
plot(dat12.res$store.param[[1]], type='l')
title(y="# of Breakpoints")

#LML
plot(dat12.res$store.param[[2]], type='l')
title(y="Log Marginal Likelihood")




### ID 19

dat19.res<- gibbs.time.seg(dat.list$`19`, 10000)
###Takes ~ 0.89 hrs to run for 10000 iterations; identified 35 breakpoints
##LML of final model is -2604.251


length(dat19.res$breakpt)
#write.csv(dat19.res$breakpt, "ID19 Breakpoints (5 km).csv", row.names = F)

image(as.matrix(obs19[,-379])); title(x="Observations (n=3068)", y="Locations (n=378)")
abline(v=breakpt/nrow(obs19),lty=3,col='blue')


## traceplots

#breakpts
plot(dat19.res$store.param[[1]], type='l')
title(y="# of Breakpoints")

#LML
plot(dat19.res$store.param[[2]], type='l')
title(y="Log Marginal Likelihood")




### ID 27

dat27.res<- gibbs.time.seg(dat.list$`27`, 10000)
###Takes ~ 0.28 hrs to run for 10000 iterations; identified 7 breakpoints
##LML of final model is -984.5128


length(dat27.res$breakpt)
#write.csv(dat27.res$breakpt, "ID27 Breakpoints (5 km).csv", row.names = F)

image(as.matrix(obs27[,-379])); title(x="Observations (n=1182)", y="Locations (n=378)")
abline(v=breakpt/nrow(obs27),lty=3,col='blue')


## traceplots

#breakpts
plot(dat27.res$store.param[[1]], type='l')
title(y="# of Breakpoints")

#LML
plot(dat27.res$store.param[[2]], type='l')
title(y="Log Marginal Likelihood")



############################
#### Import Breakpoints ####
############################

obs<- read.csv("Occupancy Matrix for all Obs and Locs.csv", header = T, sep = ",")
obs1.breakpts<- read.csv("ID1 Breakpoints (5 km).csv", header = T, sep = ",")
obs1.breakpts=obs1.breakpts[,1]
obs12.breakpts<- read.csv("ID12 Breakpoints (5 km).csv", header = T, sep = ",")
obs12.breakpts=obs12.breakpts[,1]
obs19.breakpts<- read.csv("ID19 Breakpoints (5 km).csv", header = T, sep = ",")
obs19.breakpts=obs19.breakpts[,1]
obs27.breakpts<- read.csv("ID27 Breakpoints (5 km).csv", header = T, sep = ",")
obs27.breakpts=obs27.breakpts[,1]


obs1=data.frame(obs) %>% filter(id == 1) %>% dplyr::select(-id) %>% mutate(time1=1:nrow(.))
obs12=data.frame(obs) %>% filter(id == 12) %>% dplyr::select(-id) %>% mutate(time1=1:nrow(.))
obs19=data.frame(obs) %>% filter(id == 19) %>% dplyr::select(-id) %>% mutate(time1=1:nrow(.))
obs27=data.frame(obs) %>% filter(id == 27) %>% dplyr::select(-id) %>% mutate(time1=1:nrow(.))




######################################
#### Summarize Results From Model ####
######################################

nloc=ncol(obs)-1 #remove time1

setwd("~/Documents/Snail Kite Project/Data/R Scripts/activity_center1")

obs1.seg=get.summary.stats(obs1.breakpts,obs1,nloc)
#write.csv(obs1.seg, "ID1 Seg x Loc.csv", row.names = F)
obs12.seg=get.summary.stats(obs12.breakpts,obs12,nloc)
#write.csv(obs12.seg, "ID12 Seg x Loc.csv", row.names = F)
obs19.seg=get.summary.stats(obs19.breakpts,obs19,nloc)
#write.csv(obs19.seg, "ID19 Seg x Loc.csv", row.names = F)
obs27.seg=get.summary.stats(obs27.breakpts,obs27,nloc)
#write.csv(obs27.seg, "ID27 Seg x Loc.csv", row.names = F)


######################################
#### Assign Spatial Time Segments ####
######################################

dat1<- dat %>% filter(id==1) %>% assign.time.seg(obs1.seg, obs1.breakpts, dat1)
dat12<- dat %>% filter(id==12) %>% assign.time.seg(obs12.seg, obs12.breakpts, dat12)
dat19<- dat %>% filter(id==19) %>% assign.time.seg(obs19.seg, obs19.breakpts, dat19)
dat27<- dat %>% filter(id==27) %>% assign.time.seg(obs27.seg, obs27.breakpts, dat27)

dat<- rbind(dat1, dat12, dat19, dat27)

write.csv(dat, "Snail Kite Gridded Data.csv", row.names = F)

