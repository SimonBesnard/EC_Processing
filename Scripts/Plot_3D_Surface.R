## Script to compute statistical analysis annual carbon flux
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library(emdbook)

# 1. Plot for NEP

#1.1 NEP versus precipitation and stand age
params <- c(a0=2.378735e+02,a1=-3.295241e-03, a2=-7.878536e+02, a3=-1.518790e-01,
            p0=380.0814,p1=-0.000664)
Age<- as.matrix(seq(0:320)) 
Preci<-as.matrix(seq(from=100, to=3000, by=8))

NEP_3d<-curve3d(with(as.list(params),
             (a0*(exp(a1*Age))+a2*(exp(a3*Age)))+
               ((p0*(1-exp(p1*Preci))))),
        varnames=c("Age","Preci"),
        xlim=c(0,320),ylim=c(100,3000),
        sys3d="persp")
pdf("Latex/Figures/Surface_NEP_P_Age.eps", width = 7, height = 7) # Open a new pdf file
persp(t$x, t$y, t$z, phi = 25, theta = 35, lwd=0.2,xlab = "Stand Age", ylab = "Annual precipitation", zlab = "NEP",
      xlim=c(0,320),ylim=c(100,3000))->res
points(trans3d(NEP$Stand_Age, NEP$Annual_Preci, NEP$values, pmat = res), col="black", pch=20, cex=0.5)
dev.off() # Close the file

#1.1 NEP versus stand age and temperature
params <- c(a0=2.378735e+02,a1=-3.295241e-03, a2=-7.878536e+02, a3=-1.518790e-01,
            t0=384.0151, t1=1.315, t2=-0.119)
Age<- as.matrix(seq(0:320)) 
Temp<- as.matrix(seq(from=-5, to=30, by=0.1))

NEP_3d<-curve3d(with(as.list(params),
             (a0*(exp(a1*Age))+a2*(exp(a3*Age)))+
               ((t0/(1+exp(t1+t2*Temp))))),
        varnames=c("Age","Temp"),
        xlim=c(0,320),ylim=c(-5, 30),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Temperature", zlab = "NEP",
        phi = 25, theta = 35, lwd=0.2)
pdf("Latex/Figures/Surface_NEP_T_Age.eps", width = 7, height = 7) # Open a new pdf file
persp(NEP_3d$x, NEP_3d$y, NEP_3d$z, phi = 25, theta = 35, lwd=0.2, xlab = "Stand Age", ylab = "Temperature", zlab = "NEP",
      xlim=c(0,320),ylim=c(-5,30))->res
points(trans3d(NEP$Stand_Age, NEP$Tair, NEP$values, pmat = res), col="black", pch=20, cex=0.5)
dev.off() # Close the file

# 2. Plot for GPP

#2.1 GPP versus precipitation and stand age
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

#2.1 GPP versus stand age and temperature
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

# 3. Plot for Reco

#3.1 Reco versus precipitation and stand age
params <- c(a0=631.617082164,a1=0.154251100, a2=-0.001269377,
            p0=2640.012,p1=-0.000664)
Age<- as.matrix(seq(0:350)) 
Preci<-as.matrix(seq(from=10, to=3000, by=8))

pdf("Latex/Figures/Surface_Reco_P_Age.eps", width = 7, height = 7) # Open a new pdf file
curve3d(with(as.list(params),
             a0*(Age^a1)*(exp(a2*Age))+
               ((p0*(1-exp(p1*Preci))))),
        varnames=c("Age","Preci"),
        xlim=c(0,350),ylim=c(500,3000),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Annual precipitation", zlab = "Reco",
        phi = 25, theta = 35, lwd=0.2)
dev.off() # Close the file

#3.2 Reco versus stand age and temperature
params <- c(a0=631.617082164,a1=0.154251100, a2=-0.001269377,
            t0=2487.264, t1=1.315, t2=-0.119)
Age<- as.matrix(seq(0:350)) 
Temp<- as.matrix(seq(from=-5, to=30, by=0.1))

pdf("Latex/Figures/Surface_Reco_T_Age.eps", width = 7, height = 7) # Open a new pdf file
curve3d(with(as.list(params),
             a0*(Age^a1)*(exp(a2*Age))+
               ((t0/(1+exp(t1+t2*Temp))))),
        varnames=c("Age","Temp"),
        xlim=c(0,350),ylim=c(-5, 30),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Temperature", zlab = "Respiration",
        phi = 25, theta = 35, lwd=0.2)
dev.off() # Close the file

# 4. Plot for Ratio NEP-GPP

#4.1 Ratio NEP-GPP versus precipitation and stand age
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

#4.2 Ratio NEP-GPP versus stand age and temperature
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

# 5. Plot for Ratio GPP-Reco

#5.1 Ratio GPP-Reco versus precipitation and stand age
params <- c(a0=1.1582127, a1=-0.2312178,
            p0= 1.14234501, p1=-0.00930973)
Age<- as.matrix(seq(0:350)) 
Preci<-as.matrix(seq(from=100, to=3000, by=10))

pdf("Latex/Figures/Surface_Ratio_GPP_Reco_P_Age.eps", width = 7, height = 7) # Open a new pdf file
curve3d(with(as.list(params),
             a0*(1-exp(a1*Age))+
               ((p0*(1-exp(p1*Preci))))),
        varnames=c("Age","Preci"),
        xlim=c(0,100),ylim=c(100,3000),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Annual precipitation", zlab = "Ratio GPP-ER",
        phi = 20, theta = 35, lwd=0.2)
dev.off() # Close the file

#5.2. Ratio GPP-Reco versus stand age and temperature
params <- c(a0=1.1582127, a1=-0.2312178,
            t0=1.2423580, t1=-1.3658204, t2=-0.2106015)
Age<- as.matrix(seq(0:350)) 
Temp<- as.matrix(seq(from=-5, to=30, by=0.11))

pdf("Latex/Figures/Surface_Ratio_GPP_Reco_T_Age.eps", width = 7, height = 7) # Open a new pdf file
curve3d(with(as.list(params),
             a0*(1-exp(a1*Age))+
               ((t0/(1+exp(t1+t2*Temp))))),
        varnames=c("Age","Temp"),
        xlim=c(0,350),ylim=c(-5, 30),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Temperature", zlab = "Ratio GPP-ER", 
        phi = 25, theta = 35, lwd=0.2)
dev.off() # Close the file

# 6. Plot for Ratio NEP-GPPmax

#6.1 Ratio GPP-Reco versus precipitation and stand age
params <- c(a0=-2.143196e+00, a1=-2.588261e-06, a2=2.402659e+01,a3=-5.176585e+00, a4=2.280457e+00, a5=9.254190e-01,
            p0= 0.15361297, p1=-0.00241885)
Age<- as.matrix(seq(0:350)) 
Preci<-as.matrix(seq(from=100, to=3000, by=10))

pdf("Latex/Figures/Surface_Ratio_NEP_GPPmax_P_Age.eps", width = 7, height = 7) # Open a new pdf file
curve3d(with(as.list(params),
             (a0+a1*Age^(a2)*exp(a3*Age)+a4/(1+exp(-Age*a5)))+
               ((p0*(1-exp(p1*Preci))))),
        varnames=c("Age","Preci"),
        xlim=c(0,100),ylim=c(100,3000),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Annual precipitation", zlab = "Ratio NEP-GPPclimax",
        phi = 20, theta = 35, lwd=0.2)
dev.off() # Close the file

#6.2. Ratio NEP-GPPmax versus stand age and temperature
params <- c(a0=-2.143196e+00, a1=-2.588261e-06, a2=2.402659e+01,a3=-5.176585e+00, a4=2.280457e+00, a5=9.254190e-01,
            t0=0.1848171, t1=14.7522288, t2=-6.4277945)
Age<- as.matrix(seq(0:350)) 
Temp<- as.matrix(seq(from=-5, to=30, by=0.11))

pdf("Latex/Figures/Surface_Ratio_NEP_GPPmax_T_Age.eps", width = 7, height = 7) # Open a new pdf file
curve3d(with(as.list(params),
             (a0+a1*Age^(a2)*exp(a3*Age)+a4/(1+exp(-Age*a5)))+
               ((t0/(1+exp(t1+t2*Temp))))),
        varnames=c("Age","Temp"),
        xlim=c(0,350),ylim=c(-5, 30),
        sys3d="persp",
        xlab = "Stand Age", ylab = "Temperature", zlab = "Ratio NEP-GPPclimax", 
        phi = 25, theta = 35, lwd=0.2)
dev.off() # Close the file


