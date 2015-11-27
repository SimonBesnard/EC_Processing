## Script to model NEP vs. stand age
## Author: Simon Besnard
## 12.10.2015
###################################
## Load the necessary packages
library(ggplot2)
library(gridExtra)

# Compute GPP steady state
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")
Flux_High<-dfAll_Sites[dfAll_Sites$Study %in% c("Yes"),]
GPP<- Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
params <- c(GPP15= 1923, GPP1000=1827, a1=242, a2=0.049, k=-0.00025)
GPP$GPPmat<-with(as.list(params), GPP15*((1+exp(a1-a2*15))/(1+exp(a1-a2*GPP$Tair))))
GPP$GPPp<-with(as.list(params), GPP1000*((1-exp(-k*GPP$Annual_Preci))/(1-exp(-k*1000))))
GPP<-transform(GPP, GPPmax = pmin(GPPmat, GPPp))

# Create GPP dataframe
source("Function/C_Model.R")
nyears=300
GPP_Flux<- as.data.frame(seq(0:300))
colnames(GPP_Flux)<- "Age"
f_Age_GPP<- function (x) {1287.1921959 *(1-exp(-0.1343273*x))}
GPP_Flux$f_Age<- f_Age_GPP(GPP_Flux$Age)
GPP_Flux<-GPP_Flux$f_Age

# Set parameters values of the model
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
# Run the model
Cdyn1 <- simpleCdyn(GPP_Flux, NPPeff,kcBio , kCsoil1, kCsoil2, h, 0, Csoil1SS, Csoil2SS)

# Plot outputs
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

