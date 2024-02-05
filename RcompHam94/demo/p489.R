print(Newey.West)


print(Dickey.Fuller)


print(Phillips.Perron)


print(Wald.F.Test)


data(gnptbill, package="RcompHam94")
dataset <- window(
	cbind( i=gnptbill[,"TBILL"], y=100*log(gnptbill[,"GNP"]), tt=1:dim(gnptbill)[[1]] ),
	start=c(1947,1), end=c(1989,1) )


plot(index(gnptbill),gnptbill$TBILL, type = "l", xlab="Figure 17.2 - Nominal Interest Rate", ylab="")


case1.lms <- summary( dynlm(i ~ 0 + L(i), dataset) )
case1.DF <- Dickey.Fuller( T=length(case1.lms$residuals),
  rho=case1.lms$coefficients[["L(i)","Estimate"]],
  sigma.rho=case1.lms$coefficients[["L(i)","Std. Error"]] )
print( t(case1.lms$coefficients[, c("Estimate","Std. Error"),drop=FALSE]) )
print( case1.DF )


case2.lms <- summary(dynlm( i ~ 1 + L(i), dataset))
case2.DF <- Dickey.Fuller( T=length(case2.lms$residuals),
  rho=case2.lms$coefficients[["L(i)","Estimate"]],
  sigma.rho=case2.lms$coefficients[["L(i)","Std. Error"]] )
print( t(case2.lms$coefficients[, c("Estimate","Std. Error"),drop=FALSE]) )
print( case2.DF )


F <- Wald.F.Test( R=diag(2),
                      b=case2.lms$coefficients[,"Estimate"],
                      r=c(0,1),
                      s2=case2.lms$sigma^2,
                      XtX_1=case2.lms$cov.unscaled )
print(F)


library(urca)
args(ur.df)
tbill.1.ur.df <-ur.df(dataset[,"i"],  type ="none", lags = 0)
print(summary(tbill.1.ur.df))
tbill.2.ur.df <-ur.df(dataset[,"i"],  type ="drift", lags = 0)
print(summary(tbill.2.ur.df))


print(case1.lms$coefficients["L(i)","Estimate"])
print(attr(tbill.1.ur.df,"testreg")$coefficients["z.lag.1","Estimate"]+1)
print(case2.lms$coefficients["L(i)","Estimate"])
print(attr(tbill.2.ur.df,"testreg")$coefficients["z.lag.1","Estimate"]+1)


plot(index(gnptbill),gnptbill$GNP, type = "l", xlab="Figure 17.3 - Real GNP", ylab="")


case4.lms <- summary(dynlm( y ~ 1 + L(y) + tt, dataset ))
case4.DF <- Dickey.Fuller( T=length(case4.lms$residuals),
  rho=case4.lms$coefficients[["L(y)","Estimate"]],
  sigma.rho=case4.lms$coefficients[["L(y)","Std. Error"]] )
print( t(case4.lms$coefficients[, c("Estimate","Std. Error"),drop=FALSE]) )
print( case4.DF )
F <- Wald.F.Test( R=cbind( rep(0,2), diag(2) ),
                      b=case4.lms$coefficients[,"Estimate"],
                      r=c(1,0),
                      s2=case4.lms$sigma^2,
                      XtX_1=case4.lms$cov.unscaled )
print(F)


gnp.4.ur.df <-ur.df(dataset[,"y"],  type ="trend", lags = 0)
print(summary(gnp.4.ur.df))


case2.PP <- Phillips.Perron( T=length(case2.lms$residuals),
  rho=case2.lms$coefficients[["L(i)","Estimate"]],
  sigma.rho=case2.lms$coefficients[["L(i)","Std. Error"]],
  s=case2.lms$sigma,
  lambda.hat.sq=as.numeric(Newey.West( case2.lms$residuals %o% 1, 4 )),
  gamma0=mean(case2.lms$residuals^2) )
print( t(case2.lms$coefficients[, c("Estimate","Std. Error"),drop=FALSE]) )
print( case2.PP)
case4.PP <- Phillips.Perron( T=length(case4.lms$residuals),
  rho=case4.lms$coefficients[["L(y)","Estimate"]],
  sigma.rho=case4.lms$coefficients[["L(y)","Std. Error"]],
  s=case4.lms$sigma,
  lambda.hat.sq=as.numeric(Newey.West( case4.lms$residuals %o% 1, 4 )),
  gamma0=mean(case4.lms$residuals^2) )
print( t(case4.lms$coefficients[, c("Estimate","Std. Error"),drop=FALSE]) )
print( case4.PP)


args(ur.pp)
case2.ur.pp <- ur.pp( dataset[,"i"], type ="Z-tau", model = "constant", use.lag = 4 )
print(summary(case2.ur.pp))
case4.ur.pp <- ur.pp( dataset[,"y"], type ="Z-tau", model = "trend", use.lag = 4 )
print(summary(case4.ur.pp))


tbill.lms <- summary(dynlm( i ~ L(d(i), 1:4) + 1 + L(i), dataset))
tbill.adf <- Dickey.Fuller(
  T=length(tbill.lms$residuals),
  rho=tbill.lms$coefficients[["L(i)","Estimate"]],
  sigma.rho=tbill.lms$coefficients[["L(i)","Std. Error"]],
  zeta=tbill.lms$coefficients[paste("L(d(i), 1:4)", 1:4, sep = ""),"Estimate"] )
print( t(tbill.lms$coefficients[, c("Estimate","Std. Error"),drop=FALSE]) )
print( tbill.adf )


print( tbill.lms$coefficients[["L(d(i), 1:4)4","t value"]] )


gnp.lms <- summary(dynlm( y ~ L(d(y), 1:4) + 1 + L(y) + tt, dataset))
gnp.adf <- Dickey.Fuller(
  T=length(gnp.lms$residuals),
  rho=gnp.lms$coefficients[["L(y)","Estimate"]],
  sigma.rho=gnp.lms$coefficients[["L(y)","Std. Error"]],
  zeta=gnp.lms$coefficients[paste("L(d(y), 1:4)", 1:4, sep = ""),"Estimate"] )
F <- Wald.F.Test( R=cbind( rep(0,2) %o% rep(0,5), diag(2) ),
                      b=gnp.lms$coefficients[,"Estimate"],
                      r=c(1,0),
                      s2=gnp.lms$sigma^2,
                      XtX_1=gnp.lms$cov.unscaled )
print( t(gnp.lms$coefficients[, c("Estimate","Std. Error"),drop=FALSE]) )
print( gnp.adf )
print(F)


tbill.ur.df <-ur.df(dataset[,"i"],  type ="drift", lags = 4)
print(summary(tbill.ur.df))

gnp.ur.df <-ur.df(dataset[,"y"],  type ="trend", lags = 4)
print(summary(gnp.ur.df))


t.value <- (1 - gnp.lms$coefficients[["L(y)","Estimate"]]) / gnp.lms$coefficients[["L(y)","Std. Error"]]
print( t.value )
print( (1 - pt( t.value, length(gnp.lms$residuals) )) / 2 )


for ( lag in 10:1 )
{
  gnp.lm <- dynlm( formula=as.formula(paste("y ~ L(d(y), 1:",lag,") + 1 + L(y) + tt",sep="")), data=dataset )
  if ( summary(gnp.lm)$coefficients[[paste("L(d(y), 1:",lag,")",lag,sep=""),"Pr(>|t|)"]] < .05 )
    break
}
print(lag)


