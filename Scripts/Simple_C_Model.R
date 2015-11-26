## Script to compute simple carbon model
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library (ggplot2)

# 1.Import dataframe and subset data set
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")
Flux_High<-dfAll_Sites[dfAll_Sites$Study %in% c("Yes"),]
NEP<- Flux_High[Flux_High$Type_Flux %in% c("NEP"),]
GPP<- Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
Reco<- Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]
Ratio_NEP_GPP<- Flux_High[Flux_High$Type_Flux %in% c("NEP_GPP"),]
Ratio_GPP_Reco<- Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]

# 1. Compute simple model for NEP

# Compute the piecewise linear function
y<-2.378735e+02*(exp(-3.295241e-03*300)) -7.878536e+02*(exp(-1.518790e-01*300))

# Create dataframe of the two functions
fun1 <- function(x) 60.762126 * x - 549.98
dat1 <- data.frame(x = c(0, 12.547923349), y = NA)
dat1$y <- fun1(dat1$x)
dat1$Type_Flux <- "NEP"

fun2 <- function(x) -0.431181333 * x + 217.8689299
dat2 <- data.frame(x = c(12.547923349, 320), y = NA)
dat2$y <- fun2(dat2$x)
dat2$Type_Flux <- "NEP"

# Plot functions
gg1<-ggplot(mapping = aes(x, y)) +
  geom_line(data = dat1) +
  geom_line(data = dat2)+
  ylim(-600,600)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Stand age")+
  ylab("Annual carbon flux (g.m-2.y-1)")+
  geom_vline(xintercept=12.547923349, linetype="dotted", colour="red")+
  facet_wrap(~Type_Flux)

# 2. Compute simple model for GPP

# Compute stand age function
y<-1287.1921959*(1-exp(-0.1343273*300))

# Create dataframe of the two functions
fun1 <- function(x) 59.9757*x
dat1 <- data.frame(x = c(0, 21.461300986), y = NA)
dat1$y <- fun1(dat1$x)
dat1$Type_Flux <- "GPP"

fun2 <- function(x) 0.000127273*x+1287.1538181
dat2 <- data.frame(x = c(21.461300986, 320), y = NA)
dat2$y <- fun2(dat2$x)
dat2$Type_Flux <- "GPP"

# Plot functions
gg2<-ggplot(mapping = aes(x, y)) +
  geom_line(data = dat1) +
  geom_line(data = dat2)+
  ylim(0,2000)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Stand age")+
  ylab("Annual carbon flux (g.m-2.y-1)")+
  geom_vline(xintercept=21.461300986, linetype="dotted", colour="red")+
  facet_wrap(~Type_Flux)

# 3. Compute simple model for Reco

# Compute stand age function
y<-631.617082164*(300^0.154251100)*(exp(-0.001269377*300))

# Create dataframe of the two functions
fun1 <- function(x) 48.874635*x
dat1 <- data.frame(x = c(0, 22.274238621), y = NA)
dat1$y <- fun1(dat1$x)
dat1$Type_Flux <- "Respiration"

fun2 <- function(x) -0.174*x+1092.521
dat2 <- data.frame(x = c(22.274238621, 320), y = NA)
dat2$y <- fun2(dat2$x)
dat2$Type_Flux <- "Respiration"

# Plot functions
gg3<-ggplot(mapping = aes(x, y)) +
  geom_line(data = dat1) +
  geom_line(data = dat2)+
  ylim(0,2000)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Stand age")+
  ylab("Annual carbon flux (g.m-2.y-1)")+
  geom_vline(xintercept=22.274238621, linetype="dotted", colour="red")+
  facet_wrap(~Type_Flux)

# 4. Compute simple model for Ratio GPP-Reco

# Compute stand age function
y<-1.1582127*(1-exp(-0.2312178*300))

# Create dataframe of the two functions
fun1 <- function(x) 0.074807333*x
dat1 <- data.frame(x = c(0, 15.482445835), y = NA)
dat1$y <- fun1(dat1$x)
dat1$Type_Flux <- "Ratio GPP-ER"

fun2 <- function(x) 0.000000044*x+1.1581998
dat2 <- data.frame(x = c(15.482445835, 320), y = NA)
dat2$y <- fun2(dat2$x)
dat2$Type_Flux <- "Ratio GPP-ER"

# Plot functions
gg4<-ggplot(mapping = aes(x, y)) +
  geom_line(data = dat1) +
  geom_line(data = dat2)+
  ylim(0,1.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Stand age")+
  ylab("Annual carbon flux (g.m-2.y-1)")+
  geom_vline(xintercept=15.482445835, linetype="dotted", colour="red")+
  facet_wrap(~Type_Flux)

# 6. Compute simple model for Ratio NEP-GPP

# Compute stand age function
y<-0.165451439*(exp(-0.003771699*300)) -1.319022091*(exp(-0.148502307*300))

# Create dataframe of the two functions
fun1 <- function(x) 0.091073063*x -1.153571
dat1 <- data.frame(x = c(0, 14.292213699), y = NA)
dat1$y <- fun1(dat1$x)
dat1$Type_Flux <- "Ratio NEP-GPP"

fun2 <- function(x) -0.000331454*x+0.15280189
dat2 <- data.frame(x = c(14.292213699, 320), y = NA)
dat2$y <- fun2(dat2$x)
dat2$Type_Flux <- "Ratio NEP-GPP"

# Plot functions
gg5<-ggplot(mapping = aes(x, y)) +
  geom_line(data = dat1) +
  geom_line(data = dat2)+
  ylim(-1.5, 1)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Stand age")+
  ylab("Annual carbon flux (g.m-2.y-1)")+
  geom_vline(xintercept=14.292213699, linetype="dotted", colour="red")+
  facet_wrap(~Type_Flux)

# 6. Compute simple model for Ratio NEP-GPPmax

# Compute stand age function
y<- -0.673035335*(exp(-0.089879738*50)) + 0.250875785*(exp(-0.005593567*50))

# Create dataframe of the two functions
fun1 <- function(x) 0.031870127*x -0.4221595
dat1 <- data.frame(x = c(0, 19.479811487), y = NA)
dat1$y <- fun1(dat1$x)
dat1$Type_Flux <- "Ratio NEP-GPPclimax"

fun2 <- function(x) -0.0005412*x+0.20920704
dat2 <- data.frame(x = c(19.479811487, 320), y = NA)
dat2$y <- fun2(dat2$x)
dat2$Type_Flux <- "Ratio NEP-GPPclimax"

# Plot functions
gg6<-ggplot(mapping = aes(x, y)) +
  geom_line(data = dat1) +
  geom_line(data = dat2)+
  ylim(-0.5, 0.5)+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Stand age")+
  ylab("Annual carbon flux (g.m-2.y-1)")+
  geom_vline(xintercept=19.479811487, linetype="dotted", colour="red")+
  facet_wrap(~Type_Flux)

#Plot all plots together
pdf("Latex/Figures/C_Flux_Seg.eps", width = 15, height = 10) # Open a new pdf file
grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, nrow=2) # Write the grid.arrange in the file
dev.off() # Close the file
