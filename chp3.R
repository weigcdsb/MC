############################################
####### chp3: Monte Carlo Integration ######
############################################

############################################
#### EX 3.4: first MCI
set.seed(1)

nsim<-10000
u<-runif(nsim)
den<-1:nsim


## The function to be integrated
mci.ex <- function(x){(cos(50*x)+sin(20*x))^2}
par(mfrow=c(1,3))
plot(mci.ex,
     xlim=c(0,1),ylim=c(0,4),xlab="Function",ylab="")

## The monte carlo sum
hint <- mci.ex(u)
hplot <- cumsum(hint)/den
stdh <- sqrt( cumsum(hint^2)/den - (hplot)^2)

par(new=F)
hist(hint,xlab="Generated Values of Function",
     freq=F,col="grey",breaks=30,ylab="",main="")
par(new=F)
plot(hplot,type="l",col="black",
     xlab="Mean and Standard Errors",ylab="",xlim=c(1,nsim),ylim=c(.9,1.1))
par(new=T)
plot(hplot+stdh/sqrt(den),type="l",
     col="grey",xlab="",ylab="",xlim=c(1,nsim),ylim=c(.9,1.1),lty=2)
par(new=T)
plot(hplot-stdh/sqrt(den),type="l",
     col="grey",xlab="",ylab="",xlim=c(1,nsim),ylim=c(.9,1.1),lty=2)


############################################
#### EX 3.13: Student's t-distribution
set.seed(2)

## integrands
h1 <- function(x){sqrt(abs(x/(1-x)))}
h2 <- function(x){(x>2.1)*x^5}
h3 <- function(x){(x>0)*x^5/(1+(x-3)^2)}

# Ranges of estimators for various importance functions

stu <- matrix(rt(500*2000,df=12),ncol=500)
nore <- matrix(rnorm(500*2000,sd=sqrt(1.2)),ncol=500)
cauy <- matrix(rcauchy(500*2000),ncol=500)
tgam <- matrix((1+(1-2*(runif(500*2000)<.5))*rgamma(500*2000,.5)),ncol=500)

# h1
f <- h1
so1 <- apply(apply(f(stu),2,cumsum)/(1:2000),1,range)
no1 <- apply(apply(f(nore)*dt(nore,df=12)/dnorm(nore,sd=sqrt(1.2)),2,cumsum)/(1:2000),1,range)
co1 <- apply(apply(f(cauy)*dt(cauy,df=12)/dcauchy(cauy),2,cumsum)/(1:2000),1,range)
go1 <- apply(apply(f(tgam)*dt(tgam,df=12)*2/dgamma(abs(tgam-1),.5),2,cumsum)/(1:2000),1,range)

## Plots
par(mfrow=c(1,4),mar=c(2,2,1,1))
# Student
temp <- apply(apply(f(stu),2,cumsum)/1:2000,1,mean)
plot(temp,lwd=2,col="sienna3",xlab="iterations",ylab="",type="l",ylim=c(0,3))
polygon(c(1:2000,2000:1),c(so1[1,],so1[2,2000:1]),col="gray91")
lines(temp,lwd=2,col="sienna3")

# Normal
temp <- apply(apply(f(nore)*dt(nore,df=12)/dnorm(nore,sd=sqrt(1.2)),2,cumsum)/1:2000,1,mean)
plot(temp,lwd=2,col="sienna3",xlab="iterations",ylab="",type="l",ylim=c(0,3))
polygon(c(1:2000,2000:1),c(no1[1,],no1[2,2000:1]),col="gray91")
lines(temp,lwd=2,col="sienna3")

# Cauchy
temp <- apply(apply(f(cauy)*dt(cauy,df=12)/dcauchy(cauy),2,cumsum)/1:2000,1,mean)
plot(temp,lwd=2,col="sienna3",xlab="iterations",ylab="",type="l",ylim=c(0,3))
polygon(c(1:2000,2000:1),c(co1[1,],co1[2,2000:1]),col="gray91")
lines(temp,lwd=2,col="sienna3")


# Double gamma
temp <- apply(apply(2*f(tgam)*dt(tgam,df=12)/dgamma(abs(tgam-1),.5),2,cumsum)/1:2000,1,mean)
plot(temp,lwd=2,col="sienna3",xlab="iterations",ylab="",type="l",ylim=c(0,3))
polygon(c(1:2000,2000:1),c(go1[1,],go1[2,2000:1]),col="gray91")
lines(temp,lwd=2,col="sienna3")











