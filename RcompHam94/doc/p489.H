\subsection{Preamble}
This section uses a few utility functions that follow procedures in the text for testing hypotheses about unit roots.
First is the Newey West estimator described by [10.5.10] and [10.5.15].
\begin{Scode}{}
print(Newey.West)
\end{Scode}
Next are the Dickey Fuller stats described in [17.4.7] and [17.4.9], with an optional correction for
serial correlation defined in [17.7.35] and [17.7.38].
\begin{Scode}{}
print(Dickey.Fuller)
\end{Scode}
The Phillips Perron stats are defined by [17.6.8] and [17.6.12]
\begin{Scode}{}
print(Phillips.Perron)
\end{Scode}
Finally the Wald form of an F test as defined by [8.1.32].
\begin{Scode}{}
print(Wald.F.Test)
\end{Scode}
For the following analyses we will use the R package dynlm which extends
the formula language of the workhorse lm function of R to include constructs
for expressing lags and differences.  The raw data used is a series of
treasury bill rates and real GNP.  The GNP numbers are converted to logs
and multiplied by 100 to get percentage growth rates, and we will use
data from 1947:Q1 to 1989:Q1.  Note that the text specifies a start date
of 1947:Q2, but we include Q1 because it will be used in the lag calculation
for the first "official" data point of Q2.
\begin{Scode}{}
data(gnptbill, package="RcompHam94")
dataset <- window(
	cbind( i=gnptbill[,"TBILL"], y=100*log(gnptbill[,"GNP"]), tt=1:dim(gnptbill)[[1]] ),
	start=c(1947,1), end=c(1989,1) )
\end{Scode}
\subsection{Dickey Fuller Tests for Unit Roots}
Page 489 describes the analysis of nominal
three month U.S. Treasury
yield data from
dataset gnptbill, shown below.
\begin{center}
\begin{Scode}{fig=TRUE, echo=FALSE}
plot(index(gnptbill),gnptbill$TBILL, type = "l", xlab="Figure 17.2 - Nominal Interest Rate", ylab="")
\end{Scode}
\end{center}
The regression model is shown in [17.4.13], and the results are shown below.
\begin{Scode}{}
case1.lms <- summary( dynlm(i ~ 0 + L(i), dataset) )
case1.DF <- Dickey.Fuller( T=length(case1.lms$residuals),
  rho=case1.lms$coefficients[["L(i)","Estimate"]],
  sigma.rho=case1.lms$coefficients[["L(i)","Std. Error"]] )
print( t(case1.lms$coefficients[, c("Estimate","Std. Error"),drop=FALSE]) )
print( case1.DF )
\end{Scode}
A similar analysis is described on page 494 , but a constant is included in the regression model [17.4.37].
\begin{Scode}{}
case2.lms <- summary(dynlm( i ~ 1 + L(i), dataset))
case2.DF <- Dickey.Fuller( T=length(case2.lms$residuals),
  rho=case2.lms$coefficients[["L(i)","Estimate"]],
  sigma.rho=case2.lms$coefficients[["L(i)","Std. Error"]] )
print( t(case2.lms$coefficients[, c("Estimate","Std. Error"),drop=FALSE]) )
print( case2.DF )
\end{Scode}
Example 17.5 describes how to test the joint hypothesis that the trend coefficient is 0 and the autoregressive
coefficient is 1.
\begin{Scode}{}
F <- Wald.F.Test( R=diag(2),
                      b=case2.lms$coefficients[,"Estimate"],
                      r=c(0,1),
                      s2=case2.lms$sigma^2,
                      XtX_1=case2.lms$cov.unscaled )
print(F)
\end{Scode}
We can conduct a similar analysis for cases 1 and 2 using contributed package "urca" from CRAN.
\begin{Scode}
library(urca)
args(ur.df)
tbill.1.ur.df <-ur.df(dataset[,"i"],  type ="none", lags = 0)
print(summary(tbill.1.ur.df))
tbill.2.ur.df <-ur.df(dataset[,"i"],  type ="drift", lags = 0)
print(summary(tbill.2.ur.df))
\end{Scode}
The test statistic quoted here is the t-statistic from text, and,
conveniently the critical values from the appropriate
tables in the book are printed as well.
The ur.df function uses a slightly different form of the regression than the Hamilton text, using the first difference
of the input variable as the left hand side, rather than its level.  The resulting coefficient on the lagged value
of the input variable will thus be 1 less than that obtained using the "manual" procedure above.
\begin{Scode}
print(case1.lms$coefficients["L(i)","Estimate"])
print(attr(tbill.1.ur.df,"testreg")$coefficients["z.lag.1","Estimate"]+1)
print(case2.lms$coefficients["L(i)","Estimate"])
print(attr(tbill.2.ur.df,"testreg")$coefficients["z.lag.1","Estimate"]+1)
\end{Scode}
\subsection{Analyzing GNP data}
A similar analysis can be conducted on log real GNP data described beginning on page 501, shown below.
\begin{Scode}{fig=TRUE, echo=FALSE}
plot(index(gnptbill),gnptbill$GNP, type = "l", xlab="Figure 17.3 - Real GNP", ylab="")
\end{Scode}
The regression model here incorporates a time trend, based on the shape of the GDP graph
\begin{Scode}{}
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
\end{Scode}
Similarly we can use ur.df to get similar results,
although ur.df gives a slightly different value for the intercept
because of the way the trend dummy is calculated.  
\begin{Scode}{}
gnp.4.ur.df <-ur.df(dataset[,"y"],  type ="trend", lags = 0)
print(summary(gnp.4.ur.df))
\end{Scode}
\subsection{Using Phillips Perron Tests}
Examples 17.6 and 17.7 reanalyze the case 2 and case 4 regressions above using the Phillips Perron tests
as shown on pages 511-513.
\begin{Scode}{}
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
\end{Scode}
We can get similar results for the PP t-stat by using the function ur.pp from the urca package.
Note that the implied regression generates the trend centered at T/2, so
that the intercept term shown is different than that calculated above,
where the trend starts at 1.
\begin{Scode}{}
args(ur.pp)
case2.ur.pp <- ur.pp( dataset[,"i"], type ="Z-tau", model = "constant", use.lag = 4 )
print(summary(case2.ur.pp))
case4.ur.pp <- ur.pp( dataset[,"y"], type ="Z-tau", model = "trend", use.lag = 4 )
print(summary(case4.ur.pp))
\end{Scode}
\subsection{Augmented Dickey Fuller Tests}
Example 17.8 illustrates incorporates the use of lagged regressors to (putatively) eliminate serial
correlation in the residuals.  
\begin{Scode}{}
tbill.lms <- summary(dynlm( i ~ L(d(i), 1:4) + 1 + L(i), dataset))
tbill.adf <- Dickey.Fuller(
  T=length(tbill.lms$residuals),
  rho=tbill.lms$coefficients[["L(i)","Estimate"]],
  sigma.rho=tbill.lms$coefficients[["L(i)","Std. Error"]],
  zeta=tbill.lms$coefficients[paste("L(d(i), 1:4)", 1:4, sep = ""),"Estimate"] )
print( t(tbill.lms$coefficients[, c("Estimate","Std. Error"),drop=FALSE]) )
print( tbill.adf )
\end{Scode}
The next test checks whether or not the farthest lag is different from zero, i.e. whether or not the right number
of lags are included in the equation.
\begin{Scode}{}
print( tbill.lms$coefficients[["L(d(i), 1:4)4","t value"]] )
\end{Scode}
Example 17.9 performs a similar analysis for the GNP data. 
\begin{Scode}{}
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
\end{Scode}
These calculations can be obtained using ur.df, noting that the coefficient
on the lagged value of the variable will be 1 less than the method above.
\begin{Scode}{}
tbill.ur.df <-ur.df(dataset[,"i"],  type ="drift", lags = 4)
print(summary(tbill.ur.df))

gnp.ur.df <-ur.df(dataset[,"y"],  type ="trend", lags = 4)
print(summary(gnp.ur.df))
\end{Scode}
\subsection{Example 17.10 - Bayesian Test of Autoregressive Coefficient}
Page 532 describes a test on the autoregressive coefficient that weights prior probabilities.
\begin{Scode}{}
t.value <- (1 - gnp.lms$coefficients[["L(y)","Estimate"]]) / gnp.lms$coefficients[["L(y)","Std. Error"]]
print( t.value )
print( (1 - pt( t.value, length(gnp.lms$residuals) )) / 2 )
\end{Scode}
\subsection{Determining Lag Length}
Page 530 describes an iterative process to determine the correct lag length.  This is easily expressed
in terms of the structures used above.
\begin{Scode}{}
for ( lag in 10:1 )
{
  gnp.lm <- dynlm( formula=as.formula(paste("y ~ L(d(y), 1:",lag,") + 1 + L(y) + tt",sep="")), data=dataset )
  if ( summary(gnp.lm)$coefficients[[paste("L(d(y), 1:",lag,")",lag,sep=""),"Pr(>|t|)"]] < .05 )
    break
}
print(lag)
\end{Scode}
\subsection*{Annex: R Facilities}
\textbf{Further features in R}


\fbox{
\begin{minipage}{\textwidth}
Since the tests of Dickey-Fuller and of Philips Perron, the issue of unit root testing has seen tremendous research, with hundreds of papers on the topic. For a survey of the literature, see the article of \citet{PhillipsXiao1998}, or the book of \citet{MaddalaKim1998}. 

Concerning further developments, \citet{ERS1996} used a so-called GLS detrending method to test for the presence of drift and trends, and obtain tests with higher power. Concerning the lag length selection, \citet{NgPerron2001} and \citet{PerronQu2007} introduce a new information criterion which enables a better selection of the lag length. Finally, \citet{KPSS} design a test where the null hypothesis is a stationary series (around a mean or a linear trend), while the alternative is the unit root. 

In an other direction, \citet{Hansen1995} show that by adding other related variables in the testing regression, one can obtain tests with much higher power. 
\end{minipage}
}


\fbox{
\begin{minipage}{\textwidth}
Package \pkg{urca}, well documented in the book of \citet{PfaffBook2008}, contains a number of other tests:
\begin{itemize}
\item The DF-GLS test:  \fun{ur.ers}
\item A LM test: \fun{ur.sp}
\item The KPSS test of stationarity: \fun{ur.kpss}
\item A test taking into account structural breaks: \fun{ur.za}
\end{itemize}

Package \pkg{CADFtest}, described in \citet{Lupi2009}, implements the Hansen covariate test, nesting the ADF test when no covariate is given. It offers also the choice of the lag according to the \citet{NgPerron2001} MAIC criterion. 
\end{minipage}
}
