############################################
###### chp5: Monte Carlo Optimization ######
############################################

############################################
#### EX 5.2: 1st Monte Carlo Maximization
set.seed(1)

par(mfrow=c(1,2))

## The function to be optimized
mci <- function(x){(cos(50*x)+sin(20*x))^2}
plot(mci, xlim=c(0,1),ylim=c(0,4),lwd=2)
optimize(mci, c(0, 1), tol = 0.0001, maximum=T) # 3.832543

## The monte carlo maximum
nsim<-5000
u<-runif(nsim)
max(mci(u))
plot(u,mci(u))
par(mfrow=c(1,1))

############################################
#### EX 5.5: simulated annealing of Ex 5.2
set.seed(2)

par(mfrow=c(1,2))

## The function to be optimized
mci <- function(x){(cos(50*x)+sin(20*x))^2}
plot(mci, xlim=c(0,1),ylim=c(0,4),lwd=2)
# optimize(mci, c(0, 1), tol = 0.0001, maximum=TRUE)


## Simulated annealing
nsim <- 2500
u <- runif(nsim)
xval <- array(0,c(nsim,1))
r <- .5
for(i in 2:nsim){
  test <- runif(1, min=max(xval[i-1]-r,0),max=min(xval[i-1]+r,1))
  delta <- mci(test)-mci(xval[i-1])
  
  ## T_i = 1/log(i)
  rho <- min(exp(delta*log(i)/1),1)
  xval[i] <- test*(u[i]<rho)+xval[i-1]*(u[i]>rho)
}

mci(xval[nsim])
plot(xval,mci(xval),type="l",lwd=2)
par(mfrow=c(1,1))

############################################
#### EX 5.14: censored data likelihood

## This shows complete data and censored data likelihood for normal
## generate data - x is complete, y is incomplete

n<-25
a<-4.5
theta0<-4

## data generated from normal(4,1)
xdata<-c(4.238536, 6.078472, 2.917910, 5.607990,
         6.117035, 3.509471, 4.082299, 4.097428,
         2.964908, 3.010554, 2.943618, 2.913119,
         5.867795, 5.561251, 4.075173, 4.554778,
         2.214987, 4.096038, 3.776454, 5.273564,
         4.061973, 3.996325, 3.530475, 3.552799, 3.978461)
ydata<- sort(xdata*(xdata<a)+a*(xdata>a))
m<- sum(xdata<a)
theta <- seq(theta0-2,theta0+2,length=100)

## complete data
cdlike<-function(x,mu.val=0,sd.val=1)prod(dnorm(x,mu.val,1))
cdlike.vec <- NULL
for (i in 1:length(theta)) cdlike.vec[i] <- cdlike(xdata,theta[i],sd=1)

## observed data
obslike<-function(x,mu.val=0,sd.val=1)prod(dnorm(ydata[1:m],mu.val,1))*
  prod(pnorm( ydata[(m+1):n],mu.val,1))
obslike.vec<-NULL
for (i in 1:length(theta)) obslike.vec[i] <- obslike(xdata,theta[i],sd=1)

plot(theta,obslike.vec,type="l",lwd=2,lty=2,ylab="Likelihood")
par(new=T)
plot(theta,cdlike.vec,type="l",xlab="",ylab="",xaxt="n",yaxt="n",lwd=2)
legend('topright', legend = c('obs.', 'comp.'),
       lwd = 2, lty = 2:1)

############################################
#### EX 5.17: EM for censored data (Ex 5.14)
set.seed(3)

par(mfrow=c(1,2))

## likelihood: complete & obs.
plot(theta,obslike.vec,type="l",lwd=2,lty=2,ylab="Likelihood")
par(new=T)
plot(theta,cdlike.vec,type="l",xlab="",ylab="",xaxt="n",yaxt="n",lwd=2)
legend('topright', legend = c('obs.', 'comp.'),
       lwd = 2, lty = 2:1)


## EM
nt<-50
m<- sum(xdata<a)
temp<-xdata*(xdata<a)
temp<-temp[temp!=0]
xbar<-mean(temp)
that<-array(xbar,dim=c(nt,1))
for (j in 2:nt)  {
  that[j] <-(m/n)*xbar+(1-m/n)*
    (that[j-1]+dnorm(a-that[j-1])/(1-pnorm(a-that[j-1])))
}


## now do MCEM, z=missing data, nz=size of MC sample
tmc<-array(xbar,dim=c(nt,1))
nz<-50

for (j in 2:nt)  {
  z<-array(a-1,dim=c(nz,1));
  for (k in 1:nz){while(z[k] <a) z[k] <- rnorm(1,mean=tmc[j-1],sd=1)}
  zbar<-mean(z)
  tmc[j] <-(m/n)*xbar+(1-m/n)*zbar
}

plot(that,type="l",xlim=c(0,nt),ylim=c(3.5,4.5),lwd=2,
     xlab="Iteration",ylab="EM Estimate",lty=2)
par(new=T)
plot(tmc,type="l",xlim=c(0,nt),ylim=c(3.5,4.5),xlab="",
     ylab="",xaxt="n",yaxt="n",lwd=2,lty=1)
legend('topright', legend = c('EM', 'MCEM'),
       lwd = 2, lty = 2:1)

############################################
#### EX 5.19: EM for mean mixtures
#### of normal distributions
set.seed(4)


# Production of the likelihood surface of a
# normal mean mixture
sampl <- rnorm(500)+(runif(500)<.3)*2.5
par(mfrow=c(1,2))
hist(sampl,nclass=50,xlab="x",ylab="",col="steelblue2",
     main="0.3*N(2.5,1)+ 0.7*N(0,1)",proba=T)
xz <- seq(min(sampl),max(sampl),(max(sampl)-min(sampl))/1000)
lines(xz,(.3*exp(-.5*(xz-2.5)^2)+.7*exp(-.5*xz^2))/sqrt(2*pi),
      col="red",lwd=2.2)

mu1 <- seq(-0.5,0.5,.008)
mu2 <- seq(2.0,3.0,.008)
mo1 <- mu1 %*% t(mu2/mu2)
mo2 <- (mu2/mu2) %*% t(mu2)
ca1 <- -0.5*mo1*mo1
ca2 <- -0.5*mo2*mo2
like <- 0*mo1
for (i in 1:500)
  like <- like+log(0.7*exp(ca1+sampl[i]*mo1)+0.3*exp(ca2+sampl[i]*mo2))
like <- like+.1*(ca1+ca2)
like <- like-min(like)
like <- exp(like)
image(mu1,mu2,like,xlab=expression(mu[1]),ylab=expression(mu[2]),
      col = terrain.colors(100))




sampl <- rnorm(500)+(runif(500)<.3)*1.8
mu1 <- seq(-1.5,3.5,.08)
mu2 <- seq(-1.5,3.5,.08)
mo1 <- mu1%*%t(mu2/mu2)
mo2 <- (mu2/mu2)%*%t(mu2)
ca1 <- -0.5*mo1*mo1
ca2 <- -0.5*mo2*mo2
like <- 0*mo1
for (i in 1:500){
  like <- like+log(0.7*exp(ca1+sampl[i]*mo1)+0.3*exp(ca2+sampl[i]*mo2))
}
like <- like+.1*(ca1+ca2)
like <- like-min(like)


par(mar=c(4,4,1,1))
xsur <- dev.cur()
image(mu1,mu2,like,xlab=expression(mu[1]),
      ylab=expression(mu[2]),col=heat.colors(250))
## Function delay, to check points
delay <- function(n){
  for (j in 1:(n*1000)) x <- log(4)^3
}
#
muz <- matrix(0,ncol=2,nrow=100)
# Em can start from any
# possible starting point
# but the cost is a dependence on this starting point!

par(mar=c(2,4,1,1))
xlik <- dev.cur()
plot(0,0,xlab="t",ylab="log-likelihood",xlim=c(0,20),
     ylim=c(0.3*max(like),max(like)),col="white")
khol <- c("steelblue2","sienna2","slateblue4","seagreen4","tomato3")
for (i in 1:5){
  gu1 <- mean(sampl)+1.5*rnorm(1)
  gu2 <- mean(sampl)+1.5*rnorm(1)
  muz[1,] <- (c(gu1,gu2))
  for (t in 2:100){
    
    # Allocation probabilities and averages
    probs <- 1/(1+.3*sqrt(exp(2*sampl*(gu2-gu1)-gu2^2+gu1^2))/.7)
    zeds <- probs # since zeds in (0,1)
    # EM step
    gu1 <- sum(zeds*sampl)/sum(zeds)
    gu2 <- sum((1-zeds)*sampl)/sum(1-zeds)
    muz[t,] <- (c(gu1,gu2))
    # Re-Ordering, if necessary!
    #  if (gu1>gu2){
    #  muz[t,]_(c(gu2,gu1))
    #  gu2_gu1
    #  gu1_muz[t,1]
    #  }
  }
  # Plot of the EM moves
  #
  t <- 1
  dev.set(xlik)
  locx <- max(1,min(length(mu1),round(length(mu1)*(muz[t,1]-min(mu1))/(max(mu1)-min(mu1)))))
  locy <- max(1,min(length(mu2),round(length(mu2)*(muz[t,2]-min(mu2))/(max(mu2)-min(mu2)))))
  points(0,like[locx,locy],col=khol[i],pch=20)
  for (t in 2:100){
    dev.set(xsur)
    points(muz[t-1,1],muz[t-1,2],col=khol[i],pch=20)
    delay(30)
    lines(c(muz[t-1,1],muz[t,1]),c(muz[t-1,2],muz[t,2]),col=khol[i])
    points(muz[t,1],muz[t,2],col=khol[i],pch=20)
    delay(30)
    dev.set(xlik)
    lonx <- max(1,min(length(mu1),round(length(mu1)*(muz[t,1]-min(mu1))/(max(mu1)-min(mu1)))))
    lony <- max(1,min(length(mu2),round(length(mu2)*(muz[t,2]-min(mu2))/(max(mu2)-min(mu2)))))
    lines(c(t-2,t-1),c(like[locx,locy],like[lonx,lony]),col=khol[i],lwd=1.5)
    points((t-1),like[lonx,lony],col=khol[i],pch=20)
    locx <- lonx
    locy <- lony
    if ((abs(muz[t-1,1]-muz[t,1])<.0001)&&(abs(muz[t-1,2]-muz[t,2])<.0001)) 
      break
  }
}
delay(50)
dev.set(xsur)
contour(mu1,mu2,like,levels=seq(min(like),max(like),10),add=T)








