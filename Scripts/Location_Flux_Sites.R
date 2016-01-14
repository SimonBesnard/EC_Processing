## Script to map the fluxnet sites location
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library(rgdal)
library(plyr)
library(ggplot2)
library(rworldmap)

#1. Map with the location of the fluxnet sites

# 1.1 Import site information
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")
Flux_High<-dfAll_Sites[dfAll_Sites$Study %in% c("Yes"),]
NEP<- Flux_High[Flux_High$Type_Flux %in% c("NEP"),]
Site_Location<-ddply(NEP, .(Site_ID),
                     summarise,
                     Lat= mean(Lat, na.rm=T),
                     Long=mean(Long, na.rm=T))

#1.2 Create background

# read land shapefile
+# World map
worldMap <- getMap()
worldMap_Ro<-spTransform(worldMap, CRS("+proj=wintri"))
world.points <- fortify(worldMap_Ro)
world.points$region <- world.points$id
world.df <- world.points[,c("long","lat","group", "region")]
colnames(world.df)<-c("Longitude", "Latitude", "group", "region")

# 1.3 Create Graticule and bounding box (longlat)

#Read shapefile
grat <- readOGR("Input/Flux_Site_Map", layer="ne_110m_graticules_15") 

#Change CRS
grat_robinson <- spTransform(grat, CRS("+proj=wintri"))

# convert Spdf to a dataframe
grat_robinson_df <- fortify(grat_robinson)

#Read shapefile
bbox <- readOGR("Input/Flux_Site_Map", layer="ne_110m_wgs84_bounding_box")

#Change CRS
bbox_robinson <- spTransform(bbox, CRS("+proj=wintri"))

# convert Spdf to a dataframe
bbox_robinson_df <- fortify(bbox_robinson)

# 1.4. Convert flux location to robinson CRS
Site <- SpatialPointsDataFrame(coords=Site_Location[,c("Long","Lat")],data=data.frame(Site_Location),proj4string=CRS(proj="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
Site_Robinson <- spTransform(Site, CRS("+proj=wintri"))
Site_Robinson_df <- as(Site_Robinson, "data.frame")

# 1.5 Plot fluxnet site

# create a blank ggplot theme
theme_opts <- list(theme(panel.grid.minor = element_blank(),
                                                     panel.grid.major = element_blank(),
                                                     panel.background = element_blank(),
                                                     plot.background = element_rect(fill=NULL),
                                                     panel.border = element_blank(),
                                                     axis.line = element_blank(),
                                                     axis.text.x = element_blank(),
                                                     axis.text.y = element_blank(),
                                                     axis.ticks = element_blank(),
                                                     axis.title.x = element_blank(),
                                                     axis.title.y = element_blank(),
                                                     plot.title = element_text(size=22)))

gg8<-ggplot(bbox_robinson_df, aes(long,lat, group=group)) + 
    geom_polygon(fill="lightsteelblue2") +
    geom_polygon(data=world.df, aes(Longitude, Latitude, group=group, colour=NULL), fill="darkseagreen") + 
    geom_path(data=grat_robinson_df, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
    geom_point(data=Site_Robinson_df, aes(x=Long.1,y=Lat.1, group= NULL),  shape=20, color="black", size=2)+
    coord_equal() + 
    theme_opts +
    scale_fill_manual(values=c("grey", "white"), guide="none") # change colors & remove legend
print(gg8)
ggsave(file="Latex/Figures/Sites_Location.eps", width=12.5, height=8.25, dpi=72)
