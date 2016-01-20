## Script to derive NEP map
## Author: Simon Besnard
## 12.10.2015
###################################
## Load the necessary packages
library(ncdf)
library (raster)
library (maps)
library(rasterVis)
library (ggplot2)

# Load forest age map
dir <- file.path(path, 'Forest_Age/global_forestAgeClasses_2011.nc')
Forest_Age<-read.nc(open.nc(dir))

# Plot forest age map 
r<- raster(dir, varname = "age")
plot(r)


