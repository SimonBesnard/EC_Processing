## Script to descrobe flux data
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library (ggplot2)
library(scales)
library (dplyr)
library (plyr)
library(tidyr)
library(gridExtra)

# 1.Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")
NEP<-readRDS("Output/NEP.rds")
NEP_Mean_Site<-readRDS("Output/NEP_Mean_Site.rds")
GPP<-readRDS("Output/GPP.rds")
GPP_Mean_Site<-readRDS("Output/GPP_Mean_Site.rds")
Ratio_NEP_GPP<-readRDS("Output/Ratio_NEP_GPP.rds")
Ratio_NEP_GPP_Mean_Site<- readRDS("Output/Ratio_NEP_GPP_Mean_Site.rds")

#2 Add GPP data to dataframe
NEP$GPP<- GPP$values
NEP_Mean_Site$GPP<- GPP_Mean_Site$values
Ratio_NEP_GPP$GPP<- GPP$values
Ratio_NEP_GPP_Mean_Site$GPP<- GPP_Mean_Site$values

# 3. Plot NEP vs PFTs
gg1<-ggplot(aes(y = values, x = Ecosystem), data = NEP) + 
  geom_boxplot(outlier.shape = 3, size=0.3, width=0.6)+
  xlab("PFT") + ylab("NEP [gC.m-2.y-1]")+
  theme_bw(base_size = 12, base_family = "Helvetica")
print(gg1)
ggsave(file="Latex/Figures/NEP_Var_PFT.eps", width=7, height=5, dpi=72)

# Compute NEP median/sd per PFT
NEP_PFT<- ddply(NEP, "Ecosystem",
      summarize,
      NEP=median(values, na.rm=T),
      NEP_sd= sd(values, na.rm=T))

# 4. Plot NEP vs climate
gg2<-ggplot(aes(y = values, x = Climate), data = NEP) + 
  geom_boxplot(outlier.shape = 3, size=0.3, width=0.6)+
  xlab("Climate") + ylab("NEP [gC.m-2.y-1]")+
  theme_bw(base_size = 12, base_family = "Helvetica")
print(gg2)
ggsave(file="Latex/Figures/NEP_Var_Clim.eps", width=7, height=5, dpi=72)

# Compute NEP median/sd per PFT
NEP_Clim<- ddply(NEP, "Climate",
                summarize,
                NEP=median(values, na.rm=T),
                NEP_sd= sd(values, na.rm=T))

# 4. Plot relationship between NEP and covariates

# GPP
lm_GPP<- lm(values~GPP, data=NEP)
summary(lm_GPP)
mp <- predict(lm_GPP, interval = 'conf')
df_GPP <- cbind(NEP, mp)
gg1<- ggplot(aes(y = values, x = GPP), data = df_GPP) + 
  geom_point()+
  xlab("GPP [gC.m-2.y-1]") + ylab("NEP [gC.m-2.y-1]")+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  annotate("text", label = "R-squared = 0.26, P<0.001", x = 3400, y = 1900, size =4) + 
  geom_smooth(method='lm',  se = FALSE, color = 'black') +
  geom_line(aes(y = upr), color = 'black', linetype = 2) +
  geom_line(aes(y = lwr), color = 'black', linetype = 2)+ 
  ylim(-700, 2200)
print(gg1)
ggsave(file="Latex/Figures/Cor_NEP_GPP.eps", width=8, height=6, dpi=72)

# MAT
lm_MAT<- lm(values~MAT, data=NEP)
summary(lm_MAT)
mp <- predict(lm_MAT, interval = 'conf')
df_MAT <- cbind(NEP, mp)
gg2<- ggplot(aes(y = values, x = MAT), data = df_MAT) + 
  geom_point()+
  xlab("MAT [°C]") + ylab("NEP [gC.m-2.y-1]")+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  annotate("text", label = "R-squared = 0.07, P<0.001", x = 24.3, y = 1900, size =4) + 
  geom_smooth(method='lm',  se = FALSE, color = 'black') +
  geom_line(aes(y = upr), color = 'black', linetype = 2) +
  geom_line(aes(y = lwr), color = 'black', linetype = 2)+ 
  ylim(-700, 2200)
print(gg2)
ggsave(file="Latex/Figures/Cor_NEP_MAT.eps", width=8, height=6, dpi=72)

# Stand Age
lm_Stand_Age<- lm(values~Stand_Age, data=NEP)
summary(lm_Stand_Age)
gg3<- ggplot(aes(y = values, x = Stand_Age), data = NEP) + 
  geom_point()+
  xlab("Age [years]") + ylab("NEP [gC.m-2.y-1]")+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  annotate("text", label = "R-squared = 0, P=0.68", x = 315, y = 1900, size =4)+ 
  ylim(-700, 2200)
print(gg3)
ggsave(file="Latex/Figures/Cor_NEP_Age.eps", width=8, height=6, dpi=72)

# SPI 
lm_SPI_CRU<- lm(values~SPI_CRU, data=NEP)
summary(lm_SPI_CRU)
mp <- predict(lm_SPI_CRU, interval = 'conf')
df_SPI_CRU <- cbind(NEP, mp)
gg4<- ggplot(aes(y = values, x = SPI_CRU), data = df_SPI_CRU) + 
  geom_point()+
  xlab("SPI") + ylab("NEP [gC.m-2.y-1]")+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  annotate("text", label = "R-squared = 0.03, P<0.001", x = 1.5, y = 1900, size =4) + 
  geom_smooth(method='lm',  se = FALSE, color = 'black') +
  geom_line(aes(y = upr), color = 'black', linetype = 2) +
  geom_line(aes(y = lwr), color = 'black', linetype = 2) + 
  ylim(-700, 2200)
print(gg4)
ggsave(file="Latex/Figures/Cor_NEP_SPI.eps", width=8, height=6, dpi=72)

# N deposition
lm_NHx<- lm(values~NHx, data=NEP_Mean_Site)
summary(lm_NHx)
mp <- predict(lm_NHx, interval = 'conf')
df_NHx <- cbind(NEP_Mean_Site, mp)
gg5<- ggplot(aes(y = values, x = NHx), data = df_NHx) + 
  geom_point()+
  xlab("N deposition [kg N.ha-1.y-1]") + ylab("NEP [gC.m-2.y-1]")+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  annotate("text", label = "R-squared = 0.07, P<0.01", x = 39, y = 1900, size =4) + 
  geom_smooth(method='lm',  se = FALSE, color = 'black') +
  geom_line(aes(y = upr), color = 'black', linetype = 2) +
  geom_line(aes(y = lwr), color = 'black', linetype = 2)+ 
  ylim(-700, 2200)
print(gg5)
ggsave(file="Latex/Figures/Cor_NEP_Ndepo.eps", width=8, height=6, dpi=72)
