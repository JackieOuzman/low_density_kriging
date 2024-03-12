#install.packages("terra")
library(terra)
library(sp)
#library(rgdal)

library(ggplot2)
library(raster)
library(sf)
library(automap)
library(gstat)
library(parallel)
library(foreach)
library(doParallel)
#library(rlist)
library(spatialEco)
library(gfcanalysis)
library(stars)
library(tidyverse)


set.seed(1)
headDir <- "//FSSA2-ADL/clw-share1/mallee_mod/Therese_Jackie/smart_farms/sites/Walpeup/penetrometer/Jackie processing"

## 1) Read in data
pen.dat.df <- as.data.frame(read.csv(paste0(headDir,'/20230807_Walpeup_Pen_jaxfor mapping.csv')))
pen.dat.df <- pen.dat.df %>% dplyr::mutate(temp_for_dupl = paste0(x, y))
pen.dat.df <- pen.dat.df %>% dplyr::distinct(temp_for_dupl, .keep_all = TRUE) #work around because I cant get st_unique to work


pen.dat.df <-
  pen.dat.df %>% rename(
    max_peak = "max..Peak..resistance.value.up.to.50cm" ,
    location = "location.in.the.profile.of.first.peak..up.to.50cm." ,
    depth_to_resis = "The.depth.when.resistance.first.exceeds.2.5MPa.to.depth.of.50cm" ,
    sum_area_to_resist = "Sum.the.area.of.the.curve.until.2.5MPa.is.first.reached..up.to.50cm.",
    sum_area_to_50 = "X0.50"
  )


pen.dat.sf <- st_as_sf(pen.dat.df,coords=c("x","y"),crs=28354)




#em.dat.sf <- st_unique(em.dat.sf)
str(pen.dat.sf)
#boundary <- st_read("//FSSA2-ADL/clw-share1/mallee_mod/Therese_Jackie/smart_farms/sites/Walpeup/Block_boundary/Norm_Walpeup_Proj.shp")[1]
boundary <- st_read("//FSSA2-ADL/clw-share1/mallee_mod/Therese_Jackie/smart_farms/sites/Walpeup/Block_boundary/Norm_Walpeup_proj_small_bound.shp")[1]


## 2) Convert to UTM to create grid size - known in m
# Find the best zone based on site's location
# best.crs <- utm_zone(st_bbox(em.dat.sf)[1],st_bbox(em.dat.sf)[2])
# prj.crs.char <- paste0("326",noquote(substr(best.crs,start = 1,stop = 2)))
# prj.crs.num <- as.numeric(prj.crs.char)

boundary.utm <- st_transform(boundary,28354)
pen.dat.sf.utm <- st_transform(pen.dat.sf,28354)

## 3) Create grid of points to krig onto
grid_size <- 2  # specify the grid size in meters
grid_points <- st_make_grid(boundary.utm, cellsize = c(grid_size, grid_size),what="centers")
grid.mask <- st_intersection(grid_points,boundary.utm)
plot(grid.mask)

### check the points for kriging are in the block

ggplot() +
  geom_sf(data = boundary.utm, color = "black", fill = NA) +
  geom_sf(data = pen.dat.sf.utm, color = "black", fill = NA) +
  theme_bw()+
  theme(legend.position = "none",
        axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())+
  labs(title = "",
       subtitle = "")

#################################################################################
## Variable 1 ###
names(pen.dat.sf.utm)

# "max_peak"  
# "location"  
# "depth_to_resis"   
# "sum_area_to_resist"
# "sum_area_to_50" 


rm(list=c("vgm", "map","map.sf","map.rast" , "map.rast.pred"))

## Undertake kriging
#1. fit the variogram to the spatial data and select the data clm
vgm <- autofitVariogram(`sum_area_to_50`~1,   # variable 
                        pen.dat.sf.utm[,],
                        model = c("Exp","Sph" ,"Gau","Ste")) #"Exp","Sph" ,"Gau","Ste"  #fix.values = c(0,NA,NA)
plot(vgm)

vgm <- autofitVariogram(`sum_area_to_50`~1,
                        pen.dat.sf.utm[,],
                        model = c("Exp","Sph" ,"Gau")) #"Exp","Sph" ,"Gau","Ste"  #fix.values = c(0,NA,NA)

plot(vgm)



#2. Run kriging
map <- krige(sum_area_to_50~1, # variable 
             pen.dat.sf.utm[,], 
             grid.mask,
             model = vgm$var_model, 
             nmax = 100, 
             debug.level=-1)


#3. Convert to sf object
map.sf <- st_as_sf(map)
plot(map.sf)
map.rast <- rast(st_rasterize(map.sf))
map.rast <- mask(map.rast,boundary.utm)
plot(map.rast)

map.rast.pred <- map.rast[[1]]   #Just get the predicted layer
plot(map.rast.pred)

#4. saving tiff specify the variable
writeRaster(map.rast.pred,paste0(headDir,'/low_density_kriging/Walpeup_pen_sum_area_to_50_2m_small.tif'),overwrite=T)
