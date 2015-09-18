## Script to map the fluxnet sites location
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages

#1. Map with the location of the fluxnet sites

# Import site information
Site_Info<-read.csv("Input/Potential_Sites.csv", header = TRUE)
Site_Info<- Site_Info[-c(28,36,37,49,50,51,53,63),]
Site_Info<-Site_Info[Site_Info$Intensity_Replacement %in% c("High"),]

#1.1 Land shapefile

# read land shapefile
wmap <- readOGR(dsn="Input/Flux_Site_Map", layer="ne_110m_land")

#Change CRS
wmap_robinson <- spTransform(wmap, CRS("+proj=wintri"))

# convert Spdf to a dataframe
wmap_robinson_df <- fortify(wmap_robinson)

# 7.2 Graticule and bounding box (longlat)

#Read shapefile
grat <- readOGR("Input/Flux_Site_Map", layer="ne_110m_graticules_15") 

#Change CRS
grat_robinson <- spTransform(grat, CRS("+proj=wintri"))

# convert Spdf to a dataframe
grat_robinson_df <- fortify(grat_robinson)

# 1.3 Bouding box

#Read shapefile
bbox <- readOGR("Input/Flux_Site_Map", layer="ne_110m_wgs84_bounding_box")

#Change CRS
bbox_robinson <- spTransform(bbox, CRS("+proj=wintri"))

# convert Spdf to a dataframe
bbox_robinson_df <- fortify(bbox_robinson)

# 1.4. Import Fluxnet information
Site <- SpatialPointsDataFrame(coords=Site_Info[,c("y","x")],data=data.frame(Site_Info),proj4string=CRS(proj="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
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
  geom_polygon(fill="white") +
  geom_polygon(data=wmap_robinson_df, aes(long,lat, group=group, fill=hole)) + 
  geom_path(data=grat_robinson_df, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
  geom_point(data=Site_Robinson_df, aes(x=y.1,y=x.1, group= NULL),  shape=20, color="red", size=3)+
  # labs(title="Fluxnet sites location") + 
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values=c("grey", "white"), guide="none") # change colors & remove legend
print(gg8)