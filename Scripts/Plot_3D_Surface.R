## Script to create surface response plot
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library(emdbook)

# Import NEP dataframe
Flux_High<-dfAll_Sites[dfAll_Sites$Study %in% c("Yes"),]
NEP<- Flux_High[Flux_High$Type_Flux %in% c("NEP"),]

# Plot for NEP vs age and GPP
params_Age <- c(a0=0.19945826,a1=-0.00517618,
            a2=-1.51419788,a3=-0.17645323)
Age<- as.matrix(seq(0:320)) 
Preci<-as.matrix(seq(from=40, to=3000, by=11))
Temp<- as.matrix(seq(from=-5, to=30, by=0.13))
GPP<- data.frame(Preci, Temp)
params_GPP <- c(GPP15= 1923, GPP1000=1827, a1=242, a2=0.049, k=-0.00025)
GPP$GPPmat<-with(as.list(params_GPP), GPP15*((1+exp(a1-a2*15))/(1+exp(a1-a2*GPP$Temp))))
GPP$GPPp<-with(as.list(params_GPP), GPP1000*((1-exp(-k*GPP$Preci))/(1-exp(-k*1000))))
GPP<-transform(GPP, GPPmax = pmin(GPPmat, GPPp))
GPPmax<- as.matrix(GPP$GPPmax)

NEP_Surface<-curve3d(with(as.list(params),
                     (a0*(exp(a1*Age))+a2*(exp(a3*Age)))*GPPmax),
                varnames=c("Age","GPPmax"),
                xlim=c(0,320),ylim=c(0,4100),
                sys3d="persp",
                xlab = "Stand Age", ylab = "GPP", zlab = "NEP",
                phi = 20, theta = 25, lwd=0.2)

pdf("Latex/Figures/Surface_NEP.eps", width = 7, height = 7) # Open a new pdf file
persp(NEP_Surface$x, NEP_Surface$y, NEP_Surface$z, phi = 25, theta = 35, lwd=0.2, xlab = "Stand Age", ylab = "GPP", zlab = "NEP",
      xlim=c(0,320),ylim=c(0,4100))->res
  points(trans3d(NEP$Stand_Age, GPP$GPPmax, NEP$values, pmat = res), col="red", pch=20, cex=0.5)
dev.off() # Close the file

