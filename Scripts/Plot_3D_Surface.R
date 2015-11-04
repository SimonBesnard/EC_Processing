## Script to compute statistical analysis annual carbon flux
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library(emdbook)

# 1. Plot for GPP

#1.1 GPP versus precipitation and stand age
params <- c(a0=1287.1921959,a1=-0.1343273,
            p0=3188.806,p1=-0.000664)
Age<- as.matrix(seq(0:350)) 
Preci<-as.matrix(seq(from=10, to=3000, by=8))

pdf("Latex/Figures/Surface_GPP_P_Age.eps", width = 7, height = 7) # Open a new pdf file
curve3d(with(as.list(params),
             a0*(1-exp(a1*Age))+
               ((p0*(1-exp(p1*Preci))))),
        varnames=c("Age","Preci"),
        xlim=c(0,350),ylim=c(500,3000),
          sys3d="persp",
        xlab = "Stand Age", ylab = "Annual precipitation", zlab = "GPP",
        phi = 25, theta = 35, lwd=0.2)
dev.off() # Close the file

#1.1 GPP versus stand age and temperature
params <- c(a0=1287.1921959,a1=-0.1343273,
            t0=2973.416, t1=1.315, t2=-0.119)
Age<- as.matrix(seq(0:350)) 
Temp<- as.matrix(seq(from=-5, to=30, by=0.1))

pdf("Latex/Figures/Surface_GPP_T_Age.eps", width = 7, height = 7) # Open a new pdf file
curve3d(with(as.list(params),
             a0*(1-exp(a1*Age))+
               ((t0/(1+exp(t1+t2*Temp))))),
        varnames=c("Age","Temp"),
        xlim=c(0,350),ylim=c(-5, 30),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Temperature", zlab = "GPP",
        phi = 25, theta = 35, lwd=0.2)
dev.off() # Close the file

# 2. Plot for Ratio NEP-GPP

#2.1 Ratio NEP-GPP versus precipitation and stand age
params <- c(a0=0.165451439, a1=-0.003771699, a2=-1.319022091, a3= -0.148502307,
            p0= 0.095920571,p1=-0.001707187)
Age<- as.matrix(seq(0:350)) 
Preci<-as.matrix(seq(from=100, to=3000, by=10))

pdf("Latex/Figures/Surface_Ratio_NEP_GPP_P_Age.eps", width = 7, height = 7) # Open a new pdf file
curve3d(with(as.list(params),
             (a0*(exp(a1*Age)) +a2*(exp(a3*Age)))+
               ((p0*(1-exp(p1*Preci))))),
        varnames=c("Age","Preci"),
        xlim=c(0,100),ylim=c(100,3000),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Annual precipitation", zlab = "Ratio NEP-GPP",
        phi = 20, theta = 35, lwd=0.2)
dev.off() # Close the file

#2.2 Ratio NEP-GPP versus stand age and temperature
params <- c(a0=0.165451439, a1=-0.003771699, a2=-1.319022091, a3= -0.148502307,
            t0=0.2069852, t1=1.1315, t2=-0.119)
Age<- as.matrix(seq(0:350)) 
Temp<- as.matrix(seq(from=-5, to=30, by=0.11))

pdf("Latex/Figures/Surface_Ratio_NEP_GPP_T_Age.eps", width = 7, height = 7) # Open a new pdf file
curve3d(with(as.list(params),
             (a0*(exp(a1*Age)) + a2*(exp(a3*Age)))+
               ((t0/(1+exp(t1+t2*Temp))))),
        varnames=c("Age","Temp"),
        xlim=c(0,350),ylim=c(-5, 30),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Temperature", zlab = "Ratio NEP-GPP", 
        phi = 25, theta = 35, lwd=0.2)
dev.off() # Close the file
