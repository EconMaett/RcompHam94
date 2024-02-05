data( ppp, package="RcompHam94" )
selection <- window( ppp, start=c(1973,2), end=c(1989,10) )
ppp.data <- cbind(
  p=100*log(selection[,"PZUNEW"]/selection[[1,"PZUNEW"]]),
  ner=-100*log(selection[,"EXRITL"]/selection[[1,"EXRITL"]]),
  pstar=100*log(selection[,"PC6IT"]/selection[[1,"PC6IT"]])
  )
y <- as.matrix(ppp.data)


lags <- 12
n <- dim(y)[[2]]
delta.y.lag <- embed(diff(y),lags)
X <- delta.y.lag[,-(1:n)]
T <- dim(X)[[1]]
lhs <- cbind( delta.y.lag[,1:n], y[2:(T+1),] )
aux.lm <- lm( lhs ~ 1 + X, list( lhs=lhs, X=X ) )
uv <- sapply(summary(aux.lm),FUN=function(x) { x$residuals })
u <- uv[,1:n]
v <- uv[,(n+1):(2*n)]


SigmaUU <- 1/T * t(u) %*% u
SigmaVV <- 1/T * t(v) %*% v
SigmaUV <- 1/T * t(u) %*% v
eigen.results <- eigen( solve(SigmaVV) %*% t(SigmaUV) %*% solve(SigmaUU) %*% SigmaUV)
lambda <- eigen.results$values
LRT <- -T*sum(log(1-lambda))
print(SigmaUU)
print(SigmaVV)
print(SigmaUV)
print(lambda)
print(T*log(1-lambda))
print(LRT)


ca.jo.results <- ca.jo(y, type = "eigen", ecdet = "none", K = 12,
spec="transitory", season = NULL, dumvar = NULL)
summary(ca.jo.results)


ahat1 <- eigen.results$vectors[,1]
ahat1.tilde <- ahat1 / sqrt( t(ahat1) %*% SigmaVV %*% ahat1 )
ahat1.normal <- ahat1 / ahat1[[1]]
print(ahat1)
print(ahat1.tilde)
print(ahat1.normal)


D = cbind( c(1, 0, 0), c(0, 0, 1) )
SigmaVV.tilde <- t(D) %*% SigmaVV %*% D
SigmaUV.tilde <- SigmaUV %*% D
eigen.results <- eigen( solve(SigmaVV.tilde) %*% t(SigmaUV.tilde) %*% solve(SigmaUU) %*% SigmaUV.tilde)
lambda.tilde <- eigen.results$values
h <- 1
LRT <- -T*sum(log(1-lambda[1:h])) + T*sum(log(1-lambda.tilde[1:h]))
ahat1.normal.tilde <- eigen.results$vectors[,1] / eigen.results$vectors[,1][[1]]
print(SigmaVV.tilde)
print(SigmaUV.tilde)
print(lambda.tilde)
print(T*log(1-lambda.tilde))
print(LRT)
print(ahat1.normal.tilde)


h <- 1
D = c(1, -1,-1) %o% 1
SigmaVV.tilde <- t(D) %*% SigmaVV %*% D
SigmaUV.tilde <- SigmaUV %*% D
eigen.results <- eigen( solve(SigmaVV.tilde) %*% t(SigmaUV.tilde) %*% solve(SigmaUU) %*% SigmaUV.tilde)
lambda.tilde <- eigen.results$values
LRT <- -T*sum(log(1-lambda[1:h])) + T*sum(log(1-lambda.tilde[1:h]))
print(SigmaVV.tilde)
print(SigmaUV.tilde)
print(lambda.tilde)
print(T*log(1-lambda.tilde))
print(LRT)


