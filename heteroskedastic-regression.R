
##R Jagadeesh 10550
# Bayesian heteroskedastic regression


rm(list=ls())

# Simulated data
# --------------
set.seed(12434)
n      = 100
p      = 2
q      = 2
beta   = c(9,0.9)
alpha  = c(1.0,1.0)
x      = cbind(1,rnorm(n))
w      = exp(x%*%alpha)
y     = x%*%beta + sqrt(w)*rnorm(n)
error = y-x%*%beta

plot(x[,2],y,xlab="x",ylab="y")


# prior
# -----
b0        = rep(0,p)
V0        = diag(100,p)
a0        = rep(0,p)
C0        = diag(100,p)
iV0       = solve(V0)
iV0b0     = iV0%*%b0
iC0       = solve(C0)
iC0a0     = iC0%*%a0

# Initial values
# --------------
ols   = lm(y~x-1)
b     = ols$coef
e     = y-x%*%b
a     = lm(log(e^2)~x-1)$coef

# MCMC scheme
# -----------
set.seed(23164)
sd    = 0.1
burn  = 1000
M     = 1000
LAG   = 20
niter = burn+LAG*M
bs    = matrix(0,niter,p)
as    = matrix(0,niter,q)
for (iter in 1:niter){
  print(iter)
  A  = exp(x%*%a/2)
  y1 = y/A
  x1 = x/matrix(A,n,p)  
  V1 = solve(iV0+t(x1)%*%x1)
  b1 = V1%*%(iV0b0+t(x1)%*%y1)
  b  = b1+t(chol(V1))%*%rnorm(p)
  e  = y-x%*%b
  a1 = rnorm(q,a,sd)
  la = sum(dnorm(e,0,exp(x%*%a1/2),log=TRUE))-sum(dnorm(e,0,exp(x%*%a/2),log=TRUE))
  if (log(runif(1))<la){
    a = a1
  }
  bs[iter,]   = t(b)
  as[iter,]   = t(a)
}
bs = bs[seq(burn+1,niter,by=LAG),]
as = as[seq(burn+1,niter,by=LAG),]

par(mfrow=c(2*p,3))
for (i in 1:p){
  ts.plot(bs[,i],ylab="",xlab="iteration",main=paste("beta",i,sep=""))
  abline(h=beta[i],col=2)
  acf(bs[,i],ylab="",main="")
  hist(bs[,i],prob=T,xlab="",ylab="",main="")
  abline(v=beta[i],col=2)
}
for (i in 1:p){
  ts.plot(as[,i],ylab="",xlab="iteration",main=paste("alpha",i,sep=""))
  abline(h=alpha[i],col=2)
  acf(as[,i],ylab="",main="")
  hist(as[,i],prob=T,xlab="",ylab="",main="")
  abline(v=alpha[i],col=2)
}
