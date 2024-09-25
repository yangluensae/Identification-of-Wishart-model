set.seed(100)
delta=3

N=5000 #number of draws from the proposal density

Sigma=matrix(nrow = 2,ncol=2)
Sigma[1,1]=1
Sigma[2,2]=2
Sigma[1,2]=Sigma[2,1]=0.6

Sigma_inv=solve(Sigma)

Y11=3
Y22=1

X=2*runif(N)-1 ####between 1 and -1


weight=(1-X^2)^((delta-3)/2)*exp(-Sigma_inv[1,2]*X*sqrt(Y11*Y22))
mass=sum(weight)
normalizedweight=weight/mass

sum(normalizedweight)

cumulative=cumsum(normalizedweight)
 

finalsample=rep(0, N)

for (i in 1:N)
{
  a=runif(1)
  location=1 ##we start with the beginning of the table, and we search for the right location
  while (cumulative[location]<a){location=location+1}
  finalsample[i]=X[location]
}

mean(finalsample)
Sigma[1,2]/sqrt(Sigma[1,1]*Sigma[2,2])
max(finalsample)
min(finalsample)
hist(finalsample, breaks=30, main="Histogram of R[1,2]",xlab="",ylab="",freq=FALSE)
kernel <- density(finalsample, from=-1, to=1)

# Kernel density plot
plot(kernel, lwd = 2, main = "kernel density estimator of R[1,2]",xlab="")

pdf("/Users/luyan/Library/Mobile Documents/com~apple~CloudDocs/latex papers/kernel1.pdf")
####now the distribution of the biggest eigenvalue

largeeigen=smalleigen=rep(0, N)
Sharpe12=rep(0, N)

mu1=0.001
mu2=0.0015

S1=mu1^2/Y11 ##squared sharpe ratio of market 1, does not depend on the covolatility
S2=mu2^2/Y22 ##squared sharpe ratio of market 2, does not depend on the covolatility

S12=rep(0,N) ##cross squared sharpe depends on the covolatility

for (n in 1:N){
covariance=matrix(nrow = 2,ncol=2)
covariance[1,1]=Y11
covariance[2,2]=Y22
covariance[1,2]=covariance[2,1]=sqrt(Y11*Y22)*finalsample[n]

largeeigen[n]=eigen(covariance)$values[1]
smalleigen[n]=eigen(covariance)$values[2]

inverse=solve(covariance)

S=inverse[1,1]*mu1^2+inverse[2,2]*mu2^2+2*inverse[1,2]*mu1*mu2

Sharpe12[n]=S-S1-S2
}

hist(largeeigen,breaks=30,main="Histogram of lambda1",xlab="",ylab="",freq=F)

kerneleigen <- density(largeeigen, from=3, to=4)

# Kernel density plot
plot(kerneleigen, main = "kernel density estimator of largest eigenvalue",xlab="")


hist(smalleigen)

hist(Sharpe12)
S1
S2

mean(Sharpe12)
hist(Sharpe12,breaks=30)
quantile(Sharpe12,probs=c(0.1,0.25,0.5,0.75,0.9))

sigma=sqrt(Y11+Y22-2*sqrt(Y11*Y22)*finalsample)

r=0
d1=0.5*sigma*sqrt(1/250)
d2=-d1

BS=pnorm(d1)-pnorm(d2)

quantile(BS,probs=c(0.1,0.25,0.5,0.75,0.9))
