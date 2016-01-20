## Script to model NEP vs. stand age
## Author: Simon Besnard
## 12.10.2015
###################################
## Load the necessary packages
library(gridExtra)
library (ggplot2)
library(scales)
library (dplyr)
library (plyr)
library(tidyr)
library(grid)

#1. Model NEP dynamics using multiplicative model

# Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")
NEP<-readRDS("Output/NEP.rds")
NEP_Mean_Site<-readRDS("Output/NEP_Mean_Site.rds")
GPP<-readRDS("Output/GPP.rds")
GPP_Mean_Site<-readRDS("Output/GPP_Mean_Site.rds")
Reco<-readRDS("Output/Reco.rds")
Reco_Mean_Site<-readRDS("Output/Reco_Mean_Site.rds")
Ratio_GPP_Reco<-readRDS("Output/Ratio_GPP_Reco.rds")
Ratio_GPP_Reco_Mean_Site<-readRDS("Output/Ratio_GPP_Reco_Mean_Site.rds")
Ratio_NEP_GPP<-readRDS("Output/Ratio_NEP_GPP.rds")
Ratio_NEP_GPP_Mean_Site<- readRDS("Output/Ratio_NEP_GPP_Mean_Site.rds")
Ratio_NEP_GPPmax<- readRDS("Output/Ratio_NEP_GPPmax.rds")
Ratio_NEP_GPPmax_Mean_Site<- readRDS("Output/Ratio_NEP_GPPmax_Mean_Site.rds")

# 1.2 Create NEP model using coursolle function

# 1 2.1 Create daframe with Age and GPPclimax per ecosystem

# Compute GGPclimax per climate
GPPmax<-ddply(GPP, "Climate",
              summarize,
              GPPmax= mean(GPPmax, na.rm=T))

# Compute df
NEP_model<- as.data.frame(0:320)
colnames(NEP_model)<- c("Age")
NEP_model$Tropical<- 3169
NEP_model$Dry<- 480
NEP_model$Temperate<- 1292
NEP_model$Continental<- 934
NEP_model<-gather(NEP_model, variable, values, -Age)
colnames(NEP_model)<-c("Age", "Climate", "GPPmax")

# 1.2.2 Compute NEP models

# From NEP to GPP ratio
cols <- c("Tropical"="#d7191c", "Temperate"="#fdae61", "Continental"="#abd9e9", "Dry"="#2c7bb6")
NEP_model$NEP<- (0.208727985*(exp(-0.004670439*NEP_model$Age)) -1.530989643*(exp(-0.176894536*NEP_model$Age)))*NEP_model$GPPmax
gg1<- ggplot(data=NEP_model)+
  geom_line(aes(x=Age, y=NEP, colour=Climate))+
  geom_point(data=NEP, aes(x=Stand_Age,y=values, group= NULL),  shape=20, color="black", size=2)+
  theme_bw(base_size = 16, base_family = "Helvetica")+
  xlab("Stand age (year)")+ ylab("NEP (g.m-2.y-1)")+
  theme(axis.text.x = element_text(hjust = 1),
        legend.key = element_blank(),
        legend.position="right", 
        legend.box="horizontal",
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  scale_colour_manual("Bioclimate",breaks = c("Tropical", "Temperate", "Continental", "Dry"),
                      values=cols, 
                      labels=c("Tropical", "Temperate", "Continental", "Dry"))+
  ylim(-700, 1000)

# save plot
ggsave("Latex/Figures/NEP_Model.eps", height=12, width=12)

#2. Compute C stocks and NEP dynamics using C pool model

# 2.1 Compute GPP steady state
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")
Flux_High<-dfAll_Sites[dfAll_Sites$Study %in% c("Yes"),]
GPP<- Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
params <- c(GPP15= 1923, GPP1000=1827, a1=242, a2=0.049, k=-0.00025)
GPP$GPPmat<-with(as.list(params), GPP15*((1+exp(a1-a2*15))/(1+exp(a1-a2*GPP$Tair))))
GPP$GPPp<-with(as.list(params), GPP1000*((1-exp(-k*GPP$Annual_Preci))/(1-exp(-k*1000))))
GPP<-transform(GPP, GPPmax = pmin(GPPmat, GPPp))

# 2.2 Create GPP dataframe
source("Function/C_Model.R")
nyears=300
GPP_Flux<- as.data.frame(seq(0:300))
colnames(GPP_Flux)<- "Age"
f_Age_GPP<- function (x) {1287.1921959 *(1-exp(-0.1343273*x))}
GPP_Flux$f_Age<- f_Age_GPP(GPP_Flux$Age)
GPP_Flux<-GPP_Flux$f_Age

# 2.3 Set parameters values of the model
NPPeff <- 0.5
kcBio <- rep(0.05, nyears)
kCsoil1 <- rep(0.2, nyears)
kCsoil2 <- rep(0.01, nyears)
h=0.3
#Assume steady state before: mini spin up or analytical (for each pool input/rate_constant)
CbioSS <- mean(GPP$GPPmax)*NPPeff/mean(kcBio) #not enitrely correct if kCBio varies
Csoil1SS <- CbioSS*mean(kcBio)/mean(kCsoil1)
Csoil2SS <- Csoil1SS*h*mean(kCsoil1)/mean(kCsoil2)
CtotSS<- CbioSS + Csoil1SS +  Csoil2SS

# 2.5 Run the model
Cdyn1 <- simpleCdyn(GPP_Flux, NPPeff,kcBio , kCsoil1, kCsoil2, h, 0, Csoil1SS, Csoil2SS)

# 2.6 Plot outputs
gg1<-ggplot(data = Cdyn1, aes(x = time, y = NEP)) +
  geom_line()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Stand Age (year)")+
  ylab("NEP")

gg2<-ggplot(data = Cdyn1) +
  geom_line(aes(x = time, y = Ctot, colour="Ctot"))+
  geom_line(aes(x = time, y = Cbio, colour="Cbio"))+
  geom_line(aes(x = time, y = Csoil1, colour="Csoil1"))+
  geom_line(aes(x = time, y = Csoil2, colour="Csoil2"))+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="bottom", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  xlab("Stand Age (year)")+
  ylab("Carbon stock")+
  scale_colour_manual("",breaks = c("Ctot", "Cbio", "Csoil1", "Csoil2"),
                      values=c("darkgreen", "blue","brown", "black"), 
                      labels=c("Total Carbon", "Above-ground", "Below-ground", "Deadwood"))+ 
  geom_hline(yintercept=11413.33, colour="darkgreen", linetype="dotted")+
  geom_hline(yintercept=31386.67, colour="black", linetype="dotted")+
  geom_hline(yintercept=2853.334, colour="blue", linetype="dotted")+
  geom_hline(yintercept=17120, colour="brown", linetype="dotted")

pdf("Latex/Figures/NEP_Model.eps", width = 8.87, height = 5.48) # Open a new pdf file
grid.arrange(gg1, gg2, nrow=2) # Write the grid.arrange in the file
dev.off() # Close the file

