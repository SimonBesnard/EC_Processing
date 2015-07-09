## Script to analyse IAV of annual carbon flux
## Author: Simon Besnard
## 9.07.2015
###################################
## Load the necessary packages
# install.packages('ggplot2')
# install.packages("scales") 
# install.packages("lubridate")
# install.packages("gridExtra")
# install.packages("REddyProc", repos="http://R-Forge.R-project.org", type="source")
# install.packages('RNetCDF')
# install.packages('minpack.lm')
# install.packages("dplyr")
# install.packages("testthat")
# install.packages('gtools')
# install.packages("FactoMineR")
# install.packages("leaps")
# install.packages("gamm4")
# install.packages("tidyr")
# install.packages("manipulate")
# install.packages("fitdistrplus")
# install.packages("evd")
# install.packages("flexmix")
# install.packages("mosaic")
# install.packages("bootstrap")

library(RNetCDF)
library (ggplot2)
library(scales)
library(lubridate)
library(gridExtra)
library(REddyProc)
library (dplyr)
library (plyr)
require(testthat)
library(gtools)
library(FactoMineR)
library (leaps)
library (RColorBrewer)
library(gamm4)
library(tidyr)
library(manipulate)
library(MASS)
library(fitdistrplus)
library(evd)
library(flexmix)
library(mosaic)
library(bootstrap)

#1. Create a df for all fluxnet sites dataframe

#Open all files
dir <- file.path(path, 'Fluxnet_Data/Data_4disturbed')
list <- list.files(dir, pattern=glob2rx('*.nc'), full.names=TRUE)

#Loop over the list of files
Fluxnet_Site <- list()
for(i in seq_along(list)) {
  Fluxnet_Site[[i]] = fLoadFluxNCIntoDataframe(VarList.V.s=c('NEE_f','GPP_f','Reco','NEE_fqcOK'),
                                               FileName.s=list[i], NcPackage.s = 'RNetCDF')
}

# Remove data before disturbance from the list of dataframe
Site_Date<-read.csv("Input/Potential_Sites.csv", header = TRUE)
Site_Date$Measure_Date<- as.Date(Site_Date$Measure_Date)
for (i in seq_along(Fluxnet_Site)){
  Fluxnet_Site[[i]]$DateTime=as.Date(Fluxnet_Site[[i]]$DateTime)
}

for (i in seq_along(Fluxnet_Site)){
  Fluxnet_Site[[i]]<-with(Fluxnet_Site[[i]], Fluxnet_Site[[i]][(Fluxnet_Site[[i]]$DateTime > Site_Date$Measure_Date[i]),])
}
