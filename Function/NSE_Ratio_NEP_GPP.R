stat_Ratio_NEP_GPP <- function(dat) {
  id<-unique(dat$Site_ID)
  Gamma<-c()
  Second_Poly<-c()
  Third_Poly<-c()
  Amiro<-c()
  Chen<-c()
  Coursolle<-c()
  Besnard<-c()
  Gomp<-c()
  for (i in id){
    fit1 <- try(nlsLM(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = dat[dat$Site_ID != i,], 
                    start = list(A =-1.1716, B = 0.8332, k = -0.4392), control = list(maxiter = 500)), silent=TRUE);
    Gamma[[i]] <- if (inherits(fit1, "nls")) sim = predict(fit1, newdata=dat[dat$Site_ID == i,]) else NA;
    fit2 <- try(nlsLM(values~A*Stand_Age^2+B*Stand_Age+C, data = dat[dat$Site_ID != i,], 
                    start = list(A=-1.873e-05, B=5.139e-03, C= -1.638e-01), control = list(maxiter = 500)), silent=TRUE);
    Second_Poly[[i]] <- if (inherits(fit2, "nls")) sim = predict(fit2, newdata=dat[dat$Site_ID == i,])else NA; 
    fit3 <- try(nlsLM(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data = dat[dat$Site_ID != i,], 
                    start = list(A= 2.758e-07, B=-1.261e-04, C= 1.550e-02, D=-3.741e-01), control = list(maxiter = 500)), silent=TRUE); 
    Third_Poly[[i]]<- if (inherits(fit3, "nls")) sim = predict(fit3, newdata=dat[dat$Site_ID == i,]) else NA; 
    fit4 <- try(nlsLM(values~A*(1-exp(k*Stand_Age)), data = dat[dat$Site_ID != i,], 
                    start = list(A=0.12651, k= -0.03053), control = list(maxiter = 500)), silent=TRUE);
    Amiro[[i]]<- if (inherits(fit4, "nls")) sim = predict(fit4, newdata=dat[dat$Site_ID == i,]) else NA;
    fit5 <- try(nlsLM(values~A*(1+((B*((Stand_Age/C)^D)-1)/(exp(Stand_Age/C)))), data = dat[dat$Site_ID != i,], 
                      start = list(A = 9.709e-04, B = 1.829e+02, C = 2.996e+01, D=2.573e+00), control = list(maxiter = 500)), silent=TRUE); 
    Chen[[i]]<- if (inherits(fit5, "nls")) sim = predict(fit5, newdata=dat[dat$Site_ID == i,]) else NA;
    fit6 <- try(nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = dat[dat$Site_ID != i,],
                      start = list(A=0.165450, B= -0.003772, C=-1.319022, D=-0.148503), control = list(maxiter = 500)), silent=TRUE); 
    Coursolle[[i]]<- if (inherits(fit6, "nls")) sim = predict(fit6, newdata=dat[dat$Site_ID == i,]) else NA;
    fit7<-try(nlsLM(values~A+B*Stand_Age^(C)*exp(D*Stand_Age)+E/(1+exp(-Stand_Age*H)), data =  dat[dat$Site_ID != i,], 
                start = list(A = -2.198508, B = -0.000103, C = 13.611050, D=-2.748321, E=2.312422, H=0.276304), control = list(maxiter = 500)), silent=TRUE);
    Besnard[[i]]<- if (inherits(fit7, "nls")) sim = predict(fit7, newdata=dat[dat$Site_ID == i,]) else NA;
    fit8<-try(nlsLM(values~k*exp(-exp(A-B*Stand_Age)), data =dat[dat$Site_ID != i,], 
                    start = list(A = 0, B = 0.2, k = 0.5), control = list(maxiter = 500)), silent=TRUE);
    Gomp[[i]]<- if (inherits(fit8, "nls")) sim = predict(fit8, newdata=dat[dat$Site_ID == i,]) else NA;
  }
  # list(Gamma, Second_Poly, Third_Poly, Amiro, Asympt)
  list(c(NSE(melt(Gamma)$value, dat$values, na.rm=TRUE), NSE(melt(Second_Poly)$value, dat$values, na.rm=TRUE), 
         NSE(melt(Third_Poly)$value, dat$values, na.rm=TRUE), NSE(melt(Amiro)$value, dat$values, na.rm=TRUE ), 
         NSE(melt(Chen)$value, dat$values, na.rm=TRUE ), NSE(melt(Coursolle)$value, dat$values, na.rm=TRUE),
         NSE(melt(Besnard)$value, dat$values, na.rm=TRUE), NSE(melt(Gomp)$value, dat$values, na.rm=TRUE)),
       c(cor(melt(Gamma)$value, dat$values)^2, cor(melt(Second_Poly)$value, dat$values)^2, 
         cor(melt(Third_Poly)$value, dat$values)^2, cor(melt(Amiro)$value, dat$values)^2, 
         cor(melt(Chen)$value, dat$values)^2, cor(melt(Coursolle)$value, dat$values)^2,
         cor(melt(Besnard)$value, dat$values)^2, cor(melt(Gomp)$value, dat$values)^2))
}