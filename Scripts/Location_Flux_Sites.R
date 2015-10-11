## Script to map the fluxnet sites location
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library(rworldmap)
library(dplyr)
library(ggplot2)

#1. Map with the location of the fluxnet sites

# Import site information
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")
Flux_High<-dfAll_Sites[dfAll_Sites$Study %in% c("Yes"),]
NEP<- Flux_High[Flux_High$Type_Flux %in% c("NEP"),]
Site_Location<-ddply(NEP, .(Site_ID),
                     summarise,
                     Lat= mean(Lat, na.rm=T),
                     Long=mean(Long, na.rm=T))

# World map
worldMap <- getMap()
world.points <- fortify(worldMap)
world.points$region <- world.points$id

world.df <- world.points[,c("long","lat","group", "region")]
colnames(world.df)<-c("Longitude", "Latitude", "group", "region")

Map_Site <- ggplot() + 
  geom_polygon(data = world.df, aes(x = Longitude, y = Latitude, group = group), fill="grey") +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45)+
  theme_bw(base_size = 16, base_family = "Helvetica")+
  geom_point(data=Site_Location, aes(x=Long, y=Lat))
ggsave(file="Latex/Figures/Sites_Location.eps", height = 10, width = 15)
