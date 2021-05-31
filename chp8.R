############################################
######### chp8: The Slice Sampler ##########
############################################

############################################
#### EX 8.1: simple slice sampler
set.seed(1)

f <- function(y)exp(-sqrt(y))/2
nsim<-50000
x<-array(0,dim=c(nsim,1))
x[1]<--log(runif(1))

for (i in 2:nsim)  {
  w<-runif(1,min=0,max=f(x[i-1]))
  x[i]<-runif(1,min=0,max=(-log(2*w))^2)
}

hist(x,main="Slice sample",xlim=c(0,60),
     ylim=c(0,.25),freq=F,col="steelblue3",breaks=250)
par(new=T)
plot(f, 0, 70, xlim=c(0,60),ylim=c(0,.25),
     xlab="",ylab="", col = 'red',
     xaxt="n",yaxt="n",lwd=3)

############################################
#### EX 8.2: truncated normal distribution
## N(-3, 1) in [0, 1]
set.seed(2)

xes <- seq(0,1,.01)
xes <- c(xes,1-xes)
yes <- c(exp(-(seq(0,1,.01)+3)^2/2),0*seq(0,1,.01))
plot(xes[1:101],yes[1:101],type="l",lwd=2,
     col="sienna3",xlab="x",ylab=expression(f(x)))
polygon(xes,yes,col="gold")

## Generation of the slice sampler sequence
obs <- matrix(0,ncol=2,nrow=100)
obs[1,1] <- .25
obs[1,2] <- 0.5*exp(-0.5*(obs[1,1]+3)^2)
points(obs[1,1],obs[1,2],col="sienna4",pch=20)
for (i in 2:100){
  obs[i,2] <- exp(-0.5*(obs[i-1,1]+3)^2)*runif(1)
  obs[i,1] <- runif(1)*min(1,-3+sqrt(-2*log(obs[i,2])))
}

##  Plot of the slice sampling steps
for (i in 2:10){
  if (i>2) points(obs[i-1,1],obs[i-1,2],col="steelblue3",pch=20)
  lines(c(obs[i-1,1],obs[i-1,1]),c(obs[i-1,2],obs[i,2]),col="steelblue3")
  points(obs[i-1,1],obs[i,2],col="steelblue3",pch=20)
  lines(c(obs[i-1,1],obs[i,1]),c(obs[i,2],obs[i,2]),col="steelblue3")
  points(obs[i,1],obs[i,2],col="steelblue3",pch=20)
}

par(mfrow=c(1,3))
losunif <- runif(2100)

## Plot 1
plot(xes[1:101],yes[1:101],type="l",lwd=2,col="sienna3",xlab="x",ylab=expression(f(x)))
polygon(xes,yes,col="gold")
obs <- matrix(0,ncol=2,nrow=100)
obs[1,1] <- .01
obs[1,2] <- .01
points(obs[1,1],obs[1,2],col="sienna4",pch=20)
counte <- 0
for (j in 1:100){
  for (i in 2:10){
    counte <- counte+1
    obs[i,2] <- exp(-0.5*(obs[i-1,1]+3)^2)*losunif[counte]
    counte <- counte+1
    obs[i,1] <- losunif[counte]*min(1,-3+sqrt(-2*log(obs[i,2])))
  }
  points(obs[i,1],obs[i,2],col="steelblue3",pch=20)   
}

## Plot 2
plot(xes[1:101],yes[1:101],type="l",
     lwd=2,col="sienna3",xlab="x",ylab=expression(f(x)))
polygon(xes,yes,col="gold")
obs <- matrix(0,ncol=2,nrow=100)
obs[1,1] <- .99
obs[1,2] <- .0001
points(obs[1,1],obs[1,2],col="sienna4",pch=20)
counte <- 0
for (j in 1:100){
  for (i in 2:10){
    counte <- counte+1
    obs[i,2] <- exp(-0.5*(obs[i-1,1]+3)^2)*losunif[counte]
    counte <- counte+1
    obs[i,1] <- losunif[counte]*min(1,-3+sqrt(-2*log(obs[i,2])))
  }
  points(obs[i,1],obs[i,2],col="steelblue3",pch=20)   
}

## Plot 3
plot(xes[1:101],yes[1:101],type="l",lwd=2,col="sienna3",xlab="x",ylab=expression(f(x)))
polygon(xes,yes,col="gold")
obs <- matrix(0,ncol=2,nrow=100)
obs[1,1] <- .25
obs[1,2] <- .0025
points(obs[1,1],obs[1,2],col="sienna4",pch=20)
counte <- 0
for (j in 1:100){
  for (i in 2:10){
    counte <- counte+1
    obs[i,2] <- exp(-0.5*(obs[i-1,1]+3)^2)*losunif[counte]
    counte <- counte+1
    obs[i,1] <- losunif[counte]*min(1,-3+sqrt(-2*log(obs[i,2])))
  }
  points(obs[i,1],obs[i,2],col="steelblue3",pch=20)   
}

############################################
#### EX 8.3: a 3D slice sampler
set.seed(3)

x <- rep(0,5000)
for (i in 2:5000){
  omega1 <- (1+sin(3*x[i-1])^2)*runif(1)-1
  omega2 <- (1+cos(5*x[i-1])^4)*runif(1)-1
  omega3 <- sqrt(-2*log( runif(1)*exp(-x[i-1]^2/2)))
  repeat{
    
    # |y| <= omega3, uniform --> y ~ U(-omega3, omega3)
    y <- -omega3+2*omega3*runif(1)
    if ((sin(3*y)^2>omega1)&&(cos(5*y)^4>omega2)) break  
  }
  x[i] <- y
}

plo <- hist(x,nclass=75,col="grey",proba=T,xlab="x",ylab="",main="")
labs <- seq(-3,3,.01)
dense <- (1+sin(3*labs)^2)*(1+cos(5*labs)^4)*exp(-labs^2/2)
dense <- dense*max(plo$density)/max(dense)
lines(labs,dense,col="sienna4",lwd=3)

############################################
#### EX 8.9: a poor slice sampler
set.seed(4)

par(mfrow=c(4,2))
for (d in c(1,5,10,50)){
  D <- 1/d
  x <- rep(0,1000)
  x[1] <- runif(1)
  for (i in 2:1000){
    x[i] <- runif(1)*(-log(runif(1)*exp(-x[i-1]^D)))^d
  }
  x <- x^D
  plot(x,type="l",xlab="t",ylab="u")
  acf(x,main=paste("dimension",d,sep=" "))
}





