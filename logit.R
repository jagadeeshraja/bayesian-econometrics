############################################################################################
#
#  Bayesian logistic regression
rm(list=ls())
prior = function(a,b){dnorm(a,m0[1],D0[1])*dnorm(b,m0[2],D0[2])}
like  = function(a,b){prod(exp((a+b*tc)*y)/(1+exp(a+b*tc)))}
fun   = function(b){prior(b)*like(b)}
post  = function(b){prior(b)*like(b)/py}

y    = c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0)
t    = c(53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81)
tbar = 70
tc   = t - tbar

# Prior hyperparameters
m0 = c(0,0)
C0 = c(10,10)
D0 = sqrt(C0)

# Plotting the data
# -----------------
pdf(file="oring.pdf",width=12,height=6)
plot(t,y,xlab="Temperature (in Fahrenheit)",ylab="O-ring failure",xlim=c(20,100),axes=FALSE,pch=16)
axis(2,at=c(0,1))
axis(1)
box()
text(70,0.95,"2")
text(67,0.05,"3")
text(70,0.05,"2")
text(76,0.05,"2")
dev.off()

# Prior density and likelihood function
# -------------------------------------
N  = 200
as = seq(-5,5,length=N)
bs = seq(-1,1,length=N)
p  = matrix(0,N,N)
l  = matrix(0,N,N)
po = matrix(0,N,N)
for (i in 1:N){
  for (j in 1:N){
    p[i,j] = prior(as[i],bs[j])
    l[i,j] = like(as[i],bs[j])
    po[i,j] = p[i,j]*l[i,j]
  }
}

pdf(file="priorlike-oring.pdf",width=10,height=6)
par(mfrow=c(1,2))
image(as,bs,p,xlab=expression(beta[0]),ylab=expression(beta[1]))
contour(as,bs,l,xlab=expression(beta[0]),ylab=expression(beta[1]),drawlabels=FALSE,add=TRUE)
title("Prior vs likelihood")
image(as,bs,l,xlab=expression(beta[0]),ylab=expression(beta[1]))
contour(as,bs,po,xlab=expression(beta[0]),ylab=expression(beta[1]),drawlabels=FALSE,add=TRUE)
title("Likelihood vs posterior")
dev.off()

# Sampling from p(beta|y) by SIR
# ------------------------------
set.seed(121464)
M   = 100000
M1  = 10000
as1 = rnorm(M,m0[1],D0[1])
bs1 = rnorm(M,m0[2],D0[2])
l1  = rep(0,M)
for (i in 1:M)
  l1[i] = like(as1[i],bs1[i])
w       = l1/sum(l1)
ind     = sample(1:M,size=M1,replace=TRUE,prob=w)
as2  = as1[ind]
bs2  = bs1[ind]

pdf(file="priorpost-oring.pdf",width=10,height=6)
par(mfrow=c(1,2))
plot(as1,bs1,xlim=range(as),ylim=range(bs),pch=16,xlab=expression(beta[0]),ylab=expression(beta[1]))
contour(as,bs,p,add=TRUE,col=2,lwd=2,drawlabels=FALSE)
title("Prior")
plot(as2,bs2,xlim=range(as),ylim=range(bs),pch=16,xlab=expression(beta[0]),ylab=expression(beta[1]))
contour(as,bs,po,add=TRUE,col=2,lwd=2,drawlabels=FALSE)
title("Posterior")
dev.off()

pdf(file="post-oring.pdf",width=10,height=6)
par(mfrow=c(1,2))
hist(as2,breaks=seq(min(as2),max(as2),length=25),main=expression(beta[0]),xlab="",prob=TRUE)
hist(bs2,breaks=seq(min(bs2),max(bs2),length=20),main=expression(beta[1]),xlab="",prob=TRUE)
dev.off()

# Posterior summaries
# -------------------
round(rbind(
c(mean(as2),median(as2),sqrt(var(as2)),quantile(as2,c(0.025,0.975))),
c(mean(bs2),median(bs2),sqrt(var(bs2)),quantile(bs2,c(0.025,0.975)))),4)

# Posterior predictive
# --------------------
ts = seq(20,100,by=0.5)-tbar
nt = length(ts)
ps = NULL
pp = NULL
for (i in 1:nt){
  A  = exp(as2+bs2*ts[i])
  pp = cbind(pp,A/(1+A))
}
quants = apply(pp,2,quantile,c(0.05,0.5,0.95))

pdf(file="pred-oring.pdf",width=9,height=6)
par(mfrow=c(1,1))
plot(ts+tbar,quants[2,],ylab="Probability of failure",
     xlab="Temperature (in Fahrenheit)",ylim=range(quants),type="l")
lines(ts+tbar,quants[1,],lty=2)
lines(ts+tbar,quants[3,],lty=2)
points(t,y,pch=16)
points(31,quants[1,23],col=2,pch=16)
points(31,quants[2,23],col=2,pch=16)
points(31,quants[3,23],col=2,pch=16)
abline(v=31,col=2)
dev.off()

quantile(pp[,23],c(0.025,0.5,0.975))






