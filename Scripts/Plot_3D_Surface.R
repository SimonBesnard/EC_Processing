## Script to compute statistical analysis annual carbon flux
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library(emdbook)

# 1. Plot for GPP

#1.1 GPP versus precipitation and stand age
params <- c(a0=1287.1921959,a1=-0.1343273,
            p0=3047.819,p1=-0.000664)
Age<- as.matrix(seq(0:320)) 
Preci<-as.matrix(seq(from=100, to=3000, by=8))

GPP_3d<-curve3d(with(as.list(params),
             a0*(1-exp(a1*Age))+
               ((p0*(1-exp(p1*Preci))))),
        varnames=c("Age","Preci"),
        xlim=c(0,320),ylim=c(100,3000),
          sys3d="persp",
        xlab = "Stand Age", ylab = "Annual precipitation", zlab = "GPP",
        phi = 25, theta = 35, lwd=0.2)
pdf("Latex/Figures/Surface_GPP_P_Age.eps", width = 7, height = 7) # Open a new pdf file
persp(GPP_3d$x, GPP_3d$y, GPP_3d$z, phi = 25, theta = 35, lwd=0.2, xlab = "Stand Age", ylab = "Annual precipitation", zlab = "GPP",
      xlim=c(0,320),ylim=c(100,3000))->res
points(trans3d(GPP$Stand_Age, GPP$Annual_Preci, GPP$values, pmat = res), col="red", pch=20, cex=0.5)
dev.off() # Close the file

#1.2 GPP versus stand age and temperature
params <- c(a0=1287.1921959,a1=-0.1343273,
            t0=2880.45, t1=1.315, t2=-0.119)
Age<- as.matrix(seq(0:320)) 
Temp<- as.matrix(seq(from=-5, to=30, by=0.1))

GPP_3d<-curve3d(with(as.list(params),
             a0*(1-exp(a1*Age))+
               ((t0/(1+exp(t1+t2*Temp))))),
        varnames=c("Age","Temp"),
        xlim=c(0,320),ylim=c(-5, 30),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Temperature", zlab = "GPP",
        phi = 25, theta = 35, lwd=0.2)
pdf("Latex/Figures/Surface_GPP_T_Age.eps", width = 7, height = 7) # Open a new pdf file
persp(GPP_3d$x, GPP_3d$y, GPP_3d$z, phi = 25, theta = 35, lwd=0.2, xlab = "Stand Age", ylab = "Temperature", zlab = "GPP",
      xlim=c(0,320),ylim=c(-5, 30))->res
points(trans3d(GPP$Stand_Age, GPP$Tair, GPP$values, pmat = res), col="red", pch=20, cex=0.5)
dev.off() # Close the file

#2. Plot for Reco

#3.1 Reco versus precipitation and stand age
params <- c(a0=631.617082164,a1=0.154251100, a2=-0.001269377,
            p0=2640.012,p1=-0.000664)
Age<- as.matrix(seq(0:320)) 
Preci<-as.matrix(seq(from=100, to=3000, by=8))

Reco_3d<- curve3d(with(as.list(params),
             a0*(Age^a1)*(exp(a2*Age))+
               ((p0*(1-exp(p1*Preci))))),
        varnames=c("Age","Preci"),
        xlim=c(0,320),ylim=c(100,3000),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Annual precipitation", zlab = "Reco",
        phi = 25, theta = 35, lwd=0.2)
pdf("Latex/Figures/Surface_Reco_P_Age.eps", width = 7, height = 7) # Open a new pdf file
persp(Reco_3d$x, Reco_3d$y, Reco_3d$z, phi = 25, theta = 35, lwd=0.2, xlab = "Stand Age", ylab = "Annual precipitation", zlab = "Respiration",
      xlim=c(0,320),ylim=c(100, 3000))->res
  points(trans3d(Reco$Stand_Age, Reco$Annual_Preci, Reco$values, pmat = res), col="red", pch=20, cex=0.5)
dev.off() # Close the file

#2.2 Reco versus stand age and temperature
params <- c(a0=631.617082164,a1=0.154251100, a2=-0.001269377,
            t0=2487.264, t1=1.315, t2=-0.119)
Age<- as.matrix(seq(0:320)) 
Temp<- as.matrix(seq(from=-5, to=30, by=0.1))

Reco_3d<-curve3d(with(as.list(params),
             a0*(Age^a1)*(exp(a2*Age))+
               ((t0/(1+exp(t1+t2*Temp))))),
        varnames=c("Age","Temp"),
        xlim=c(0,320),ylim=c(-5, 30),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Temperature", zlab = "Respiration",
        phi = 25, theta = 35, lwd=0.2)
pdf("Latex/Figures/Surface_Reco_T_Age.eps", width = 7, height = 7) # Open a new pdf file
persp(Reco_3d$x, Reco_3d$y, Reco_3d$z, phi = 25, theta = 35, lwd=0.2, xlab = "Stand Age", ylab = "Temperature", zlab = "Respiration",
      xlim=c(0,320),ylim=c(-5, 30))->res
points(trans3d(Reco$Stand_Age, Reco$Tair, Reco$values, pmat = res), col="red", pch=20, cex=0.5)
dev.off() # Close the file

#3. Plot for Ratio NEP-GPP

#3.1 Ratio NEP-GPP versus precipitation and stand age
params <- c(a0=0.165451439, a1=-0.003771699, a2=-1.319022091, a3= -0.148502307,
            p0= 0.095920571,p1=-0.001707187)
Age<- as.matrix(seq(0:320)) 
Preci<-as.matrix(seq(from=100, to=3000, by=10))

Ratio_NEP_GPP_3d<-curve3d(with(as.list(params),
             (a0*(exp(a1*Age)) +a2*(exp(a3*Age)))+
               ((p0*(1-exp(p1*Preci))))),
        varnames=c("Age","Preci"),
        xlim=c(0,320),ylim=c(100,3000),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Annual precipitation", zlab = "Ratio NEP-GPP",
        phi = 20, theta = 35, lwd=0.2)
pdf("Latex/Figures/Surface_Ratio_NEP_GPP_P_Age.eps", width = 7, height = 7) # Open a new pdf file
persp(Ratio_NEP_GPP_3d$x, Ratio_NEP_GPP_3d$y, Ratio_NEP_GPP_3d$z, phi = 25, theta = 35, lwd=0.2, xlab = "Stand Age", ylab = "Annual precipitation", zlab = "Ratio NEP-GPP",
      xlim=c(0,320),ylim=c(10, 3000))->res
points(trans3d(Ratio_NEP_GPP$Stand_Age, Ratio_NEP_GPP$Annual_Preci, Ratio_NEP_GPP$values, pmat = res), col="red", pch=20, cex=0.5)
dev.off() # Close the file

#3.2 Ratio NEP-GPP versus stand age and temperature
params <- c(a0=0.165451439, a1=-0.003771699, a2=-1.319022091, a3= -0.148502307,
            t0=0.2069852, t1=1.1315, t2=-0.119)
Age<- as.matrix(seq(0:320)) 
Temp<- as.matrix(seq(from=-5, to=30, by=0.11))

Ratio_NEP_GPP_3d<-curve3d(with(as.list(params),
             (a0*(exp(a1*Age)) + a2*(exp(a3*Age)))+
               ((t0/(1+exp(t1+t2*Temp))))),
        varnames=c("Age","Temp"),
        xlim=c(0,320),ylim=c(-5, 30),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Temperature", zlab = "Ratio NEP-GPP", 
        phi = 25, theta = 35, lwd=0.2)
pdf("Latex/Figures/Surface_Ratio_NEP_GPP_T_Age.eps", width = 7, height = 7) # Open a new pdf file
persp(Ratio_NEP_GPP_3d$x, Ratio_NEP_GPP_3d$y, Ratio_NEP_GPP_3d$z, phi = 25, theta = 35, lwd=0.2, xlab = "Stand Age", ylab = "Temperature", zlab = "Ratio NEP-GPP",
      xlim=c(0,320),ylim=c(-5, 30))->res
points(trans3d(Ratio_NEP_GPP$Stand_Age, Ratio_NEP_GPP$Tair, Ratio_NEP_GPP$values, pmat = res), col="red", pch=20, cex=0.5)
dev.off() # Close the file