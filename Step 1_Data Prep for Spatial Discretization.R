### Data Prep for Time Segmentation Model ###

setwd("~/Documents/Snail Kite Project/Data")

library(dplyr)
library(ggplot2)
library(lubridate)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(sp)
library(raster)
library(rgdal)


#####################
#### Import Data ####
#####################

dat<- read.csv("Modified Snail Kite Data.csv", header = T, sep = ",")

dat$id<- as.factor(dat$id)
dat$ESTtime<- as_datetime(dat$ESTtime)
dat$dist<- dat$dist/1000  #convert to km
names(dat)[7]<- "dist_km"
dat$R2n<- dat$R2n/1000/1000  #convert to km^2
names(dat)[9]<- "R2n_km2"


#Turn df into SpatialPointsDataFrame
dat.coords<- data.frame(dat[,c("utmlong","utmlat")])
utm.crs<- CRS("+init=epsg:32617")
dat.spdf<- SpatialPointsDataFrame(coords = dat.coords, data = dat, proj4string = utm.crs)

#Load world map data
world <- ne_countries(scale = "medium", returnclass = "sf")
usa<- world %>% filter(admin == "United States of America") %>% st_transform(proj4string(dat.spdf))

#to define extent of plot based on all individuals in SF
dat.sf<- st_as_sf(dat.spdf)



##########################################################
#### Create Grid to Discretize Space ####
##########################################################

# 5 km w 1 cell buffer
grid_5<- raster(extent(dat.spdf) + 10000)
res(grid_5)<- 5000
proj4string(grid_5)<- utm.crs

grid_5[]<- 0
dat.spdf@data$grid.cell<- cellFromXY(grid_5, dat.spdf) #378 cells occupied from possible 4272 cells

### Write to CSV for further analysis
write.csv(dat.spdf@data, "Snail Kite Gridded Data.csv", row.names = F)


### Plot all points over grid

#Create grid cell borders
borders_5<- rasterToPolygons(grid_5, dissolve = F)
borders_5f<- fortify(borders_5)

grid_5f<- as.data.frame(grid_5, xy = TRUE)

ggplot() +
  geom_sf(data = usa) +
  coord_sf(datum = NA) +
  coord_sf(xlim = c(min(dat$utmlong-20000), max(dat$utmlong+20000)),
           ylim = c(min(dat$utmlat-20000), max(dat$utmlat+20000)), expand = FALSE) +
  geom_tile(data=grid_5f, aes(x=x, y=y), fill="transparent") +
  geom_path(data = borders_5f, aes(x=long, y=lat, group=group), size=0.25) +
  geom_point(data = dat.spdf@data, aes(x=utmlong, y=utmlat, color=id), size=0.5, alpha=0.5) +
  scale_color_viridis_d(alpha = 0.6) +
  labs(x = "Easting", y = "Northing")
