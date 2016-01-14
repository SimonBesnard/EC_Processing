stat_NEP <- function(dat) {
  id<-unique(dat$Site_ID)
  Gamma<-c()
  Second_Poly<-c()
  Third_Poly<-c()
  Amiro<-c()
  Chen<-c()
  Coursolle<-c()
  Besnard<-c()
  for (i in id){
    fit1 <- try(nlsLM(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = dat[dat$Site_ID != i,], 
                    start = list(A = -786.6097, B = 0.3989, k = -0.4694), control = list(maxiter = 500)), silent=TRUE);
    Gamma[[i]] <- if (inherits(fit1, "nls")) sim = predict(fit1, newdata=dat[dat$Site_ID == i,]) else NA;
    fit2 <- try(nlsLM(values~A*Stand_Age^2+B*Stand_Age+C, data = dat[dat$Site_ID != i,], 
                    start = list(A=-0.01035, B= 2.32765, C= 77.22196), control = list(maxiter = 500)), silent=TRUE);
    Second_Poly[[i]] <- if (inherits(fit2, "nls")) sim = predict(fit2, newdata=dat[dat$Site_ID == i,])else NA; 
    fit3 <- try(nlsLM(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data = dat[dat$Site_ID != i,], 
                    start = list(A= 0.000235, B=-0.096520, C=10.214905, D=-74.366279), control = list(maxiter = 500)), silent=TRUE); 
    Third_Poly[[i]]<- if (inherits(fit3, "nls")) sim = predict(fit3, newdata=dat[dat$Site_ID == i,]) else NA; 
    fit4 <- try(nlsLM(values~A*(1-exp(k*Stand_Age)), data = dat[dat$Site_ID != i,], 
                    start = list(A=192.93829, k=-0.08976), control = list(maxiter = 500)), silent=TRUE);
    Amiro[[i]]<- if (inherits(fit4, "nls")) sim = predict(fit4, newdata=dat[dat$Site_ID == i,]) else NA;
    fit5 <- try(nlsLM(values~A*(1+((B*((Stand_Age/C)^D)-1)/(exp(Stand_Age/C)))), data = dat[dat$Site_ID != i,], 
                       start = list(A = 2.280, B = 270.887, C = 56.420, D=1.008), control = list(maxiter = 500)), silent=TRUE); 
    Chen[[i]]<- if (inherits(fit5, "nls")) sim = predict(fit5, newdata=dat[dat$Site_ID == i,]) else NA;
    fit6 <- try(nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = dat[dat$Site_ID != i,],
                      start = list(A= 2.867e+02, B=-4.941e-03, C=-1.014e+03, D=-1.802e-01), control = list(maxiter = 500)), silent=TRUE); 
    Coursolle[[i]]<- if (inherits(fit6, "nls")) sim = predict(fit6, newdata=dat[dat$Site_ID == i,]) else NA;
    fit7<-try(nlsLM(values~A+B*Stand_Age^(C)*exp(D*Stand_Age)+E/(1+exp(-Stand_Age*H)), data =  dat[dat$Site_ID != i,], 
                start = list(A = -6.294e+01, B =5.786e+02, C = -1.661e-01, D=-1.386e-03, E=-2.190e+03, H=-2.194e-01), control = list(maxiter = 500)), silent=TRUE);
    Besnard[[i]]<- if (inherits(fit7, "nls")) sim = predict(fit7, newdata=dat[dat$Site_ID == i,]) else NA;
    }
  # list(Gamma, Second_Poly, Third_Poly, Amiro, Asympt)
  list(c(NSE(melt(Gamma)$value, dat$values, na.rm=TRUE), NSE(melt(Second_Poly)$value, dat$values, na.rm=TRUE), 
         NSE(melt(Third_Poly)$value, dat$values, na.rm=TRUE), NSE(melt(Amiro)$value, dat$values, na.rm=TRUE ), 
         NSE(melt(Chen)$value, dat$values, na.rm=TRUE ), NSE(melt(Coursolle)$value, dat$values, na.rm=TRUE),
         NSE(melt(Besnard)$value, dat$values, na.rm=TRUE)),
       c(cor(melt(Gamma)$value, dat$values)^2, cor(melt(Second_Poly)$value, dat$values)^2, 
         cor(melt(Third_Poly)$value, dat$values)^2, cor(melt(Amiro)$value, dat$values)^2, 
         cor(melt(Chen)$value, dat$values)^2, cor(melt(Coursolle)$value, dat$values)^2,
         cor(melt(Besnard)$value, dat$values)^2))
}