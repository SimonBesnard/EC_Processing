stat_GPP <- function(dat) {
  id<-unique(dat$Site_ID)
  Gamma<-c()
  Second_Poly<-c()
  Third_Poly<-c()
  Amiro<-c()
  Chen<-c()
  Coursolle<-c()
  Besnard<-c()
  Gomp<-c()
  Gomp_der<-c()
  for (i in id){
    fit1 <- try(nlsLM(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = dat[dat$Site_ID != i,], 
                    start = list(A = 411.796734, B = 0.363362, k = -0.005017), control = list(maxiter = 500)), silent=TRUE);
    Gamma[[i]] <- if (inherits(fit1, "nls")) sim = predict(fit1, newdata=dat[dat$Site_ID == i,]) else NA;
    fit2 <- try(nlsLM(values~A*Stand_Age^2+B*Stand_Age, data = dat[dat$Site_ID != i,], 
                    start = list(A=-0.0849, B=23.5815), control = list(maxiter = 500)), silent=TRUE);
    Second_Poly[[i]] <- if (inherits(fit2, "nls")) sim = predict(fit2, newdata=dat[dat$Site_ID == i,])else NA; 
    fit3 <- try(nlsLM(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age, data = dat[dat$Site_ID != i,], 
                    start = list(A=-0.0008, B=-0.0001, C= 25), control = list(maxiter = 500)), silent=TRUE); 
    Third_Poly[[i]]<- if (inherits(fit3, "nls")) sim = predict(fit3, newdata=dat[dat$Site_ID == i,]) else NA; 
    fit4 <- try(nlsLM(values~A*(1-exp(k*Stand_Age)), data = dat[dat$Site_ID != i,], 
                    start = list(A=1287.1816, k= -0.1344), control = list(maxiter = 500)), silent=TRUE);
    Amiro[[i]]<- if (inherits(fit4, "nls")) sim = predict(fit4, newdata=dat[dat$Site_ID == i,]) else NA;
    fit5 <- try(nlsLM(values~A*(1+((B*((Stand_Age/C)^D)-1)/(exp(Stand_Age/C)))), data = dat[dat$Site_ID != i,], 
                      start = list(A = 1.288e+03, B = 5.793e-04, C = 7.557e+00, D=2.883e+00), control = list(maxiter = 500)), silent=TRUE); 
    Chen[[i]]<- if (inherits(fit5, "nls")) sim = predict(fit5, newdata=dat[dat$Site_ID == i,]) else NA;
    fit6 <- try(nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = dat[dat$Site_ID != i,],
                      start = list(A=1.410e+03, B= -1.035e-03, C=-1.389e+03, D=-1.107e-01), control = list(maxiter = 500)), silent=TRUE); 
    Coursolle[[i]]<- if (inherits(fit6, "nls")) sim = predict(fit6, newdata=dat[dat$Site_ID == i,]) else NA;
    fit7<-try(nlsLM(values~A+B*Stand_Age^(C)*exp(D*Stand_Age)+E/(1+exp(-Stand_Age*H)), data =  dat[dat$Site_ID != i,], 
                start = list(A = -839.33488, B = 7.02766, C =  1.21012, D=-0.02001, E=1882.40531, H=0.20122), control = list(maxiter = 500)), silent=TRUE);
    Besnard[[i]]<- if (inherits(fit7, "nls")) sim = predict(fit7, newdata=dat[dat$Site_ID == i,]) else NA;
    fit8<-try(nlsLM(values~k*exp(-exp(A-B*Stand_Age)), data =dat[dat$Site_ID != i,], 
                    start = list(A =  0.9284, B = 0.1904, k = 1282.4569), control = list(maxiter = 500)), silent=TRUE);
    Gomp[[i]]<- if (inherits(fit8, "nls")) sim = predict(fit8, newdata=dat[dat$Site_ID == i,]) else NA;
    fit9<-try(nlsLM(values~B*(k*exp(-exp(A-B*Stand_Age)))*exp(A-B*Stand_Age), data = dat[dat$Site_ID != i,],
          start = list(A = 8.771e-01, B = 9.728e-03, k = 3.862e+05), control = list(maxiter = 500)), silent = TRUE);
    Gomp_der[[i]]<- if (inherits(fit9, "nls")) sim = predict(fit9, newdata=dat[dat$Site_ID == i,]) else NA;
    
  }
  # list(Gamma, Second_Poly, Third_Poly, Amiro, Asympt)
  list(c(NSE(melt(Gamma)$value, dat$values, na.rm=TRUE), NSE(melt(Second_Poly)$value, dat$values, na.rm=TRUE), 
         NSE(melt(Third_Poly)$value, dat$values, na.rm=TRUE), NSE(melt(Amiro)$value, dat$values, na.rm=TRUE ), 
         NSE(melt(Chen)$value, dat$values, na.rm=TRUE ), NSE(melt(Coursolle)$value, dat$values, na.rm=TRUE),
         NSE(melt(Besnard)$value, dat$values, na.rm=TRUE), NSE(melt(Gomp)$value, dat$values, na.rm=TRUE), 
         NSE(melt(Gomp_der)$value, dat$values, na.rm=TRUE)),
       c(cor(melt(Gamma)$value, dat$values)^2, cor(melt(Second_Poly)$value, dat$values)^2, 
         cor(melt(Third_Poly)$value, dat$values)^2, cor(melt(Amiro)$value, dat$values)^2, 
         cor(melt(Chen)$value, dat$values)^2, cor(melt(Coursolle)$value, dat$values)^2,
         cor(melt(Besnard)$value, dat$values)^2, cor(melt(Gomp)$value, dat$values)^2, 
         cor(melt(Gomp_der)$value, dat$values)^2))
}