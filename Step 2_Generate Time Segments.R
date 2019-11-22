#Analyze Snail Kite Data with Bayesian Partitioning Model

library(tidyverse)
library(progress)
library(furrr)
library(tictoc)
library(viridis)

source('gibbs functions2.R')
source('helper functions.R')
source('gibbs sampler2.R')


###############################
#### Load and Prepare Data ####
###############################

dat<- read.csv("Snail Kite Gridded Data_large.csv", header = T, sep = ",")
dat.list<- df.to.list(dat = dat)
dat.list<- lapply(dat.list, function(x) x %>% mutate(time1 = 1:nrow(x)))  #add row for obs number

dat.long<- map_dfr(dat.list, `[`) %>% dplyr::select(id, grid.cell, time1)  #create DF
names(dat.long)[2]<- "loc.id"
dat.long$loc.id<- dat.long$loc.id %>% factor()
levels(dat.long$loc.id)<- 1:length(unique(dat.long$loc.id))  #change from raw to modified cell ID
dat.long$loc.id<- dat.long$loc.id %>% as.character() %>% as.numeric()

dat.list2<- df.to.list(dat.long)



#######################################
#### Run Gibbs Sampler for all IDs ####
#######################################

identity<- unique(dat.long$id)
ngibbs = 10000

## Run Gibbs sampler
plan(multiprocess)  #run all MCMC chains in parallel
                    #select "multiprocess" if Unix or macOS & "multisession" if Windows
                    #refer to future::plan() for more details

dat.res<- space_segment(data = dat.list2, identity = identity, ngibbs = ngibbs)
###Takes 11.41 min to run for 10000 iterations for all IDs


#Number of breakpoints by ID
dat.res$brkpts[,-1] %>% apply(1, function(x) sum(!is.na(x)))


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)

## Heatmaps
heatmap(data = dat.list2, identity = identity, dat.res = dat.res)


######################################
#### Assign Spatial Time Segments ####
######################################

dat_out<- map(dat.list, assign.time.seg) %>% map_dfr(`[`)  #assign time seg and make as DF

setwd("~/Documents/Snail Kite Project/Data/R Scripts/activity_center1")
write.csv(dat, "Snail Kite Gridded Data_large.csv", row.names = F)

