stat1 <- function(dat) {
  id<-nrow(dat)
  Gamma<-c()
  Second_Poly<-c()
  Third_Poly<-c()
  Amiro<-c()
  Asymp<-c()
  for (i in 1:id){
    fit1 <- try(nls(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = dat[-i,], 
                    start = list(A = 100, B = 0.170, k = -0.00295)), silent=TRUE);
    Gamma[i] <- if (inherits(fit1, "nls")) sim = predict(fit1, newdata=dat[i,]) else NA;
    fit2 <- try(nls(values~A*Stand_Age^2+B*Stand_Age+C, data = dat[-i,], 
                    start = list(A=-0.4, B=50, C= 300)), silent=TRUE);
    Second_Poly[i] <- if (inherits(fit2, "nls")) sim = predict(fit2, newdata=dat[i,])else NA; 
    fit3 <- try(nls(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data = dat[-i,], 
                    start = list(A=0.02, B=-0.6, C= 50, D=200)), silent=TRUE); 
    Third_Poly[i]<- if (inherits(fit3, "nls")) sim = predict(fit3, newdata=dat[i,]) else NA; 
    fit4 <- try(nls(values~A*(1-exp(k*Stand_Age)), data = dat[-i,], 
                    start = list(A=1000, k= -0.224)), silent=TRUE); 
    Amiro[i]<- if (inherits(fit4, "nls")) sim = predict(fit4, newdata=dat[i,]) else NA;
    fit5 <- try(nls(values~SSasympOff(Stand_Age, A, lrc, c0), data = dat[-i,]), silent=TRUE); 
    Asymp[i]<- if (inherits(fit5, "nls")) sim = predict(fit5, newdata=dat[i,]) else NA
  }
  # list(Gamma, Second_Poly, Third_Poly, Amiro, Asymp)
  c(NSE(Gamma, dat$values, na.rm=TRUE), NSE(Second_Poly, dat$values, na.rm=TRUE), NSE(Third_Poly, dat$values, na.rm=TRUE), 
    NSE(Amiro, dat$values, na.rm=TRUE ), NSE(Asymp, dat$values, na.rm=TRUE ) )
}
