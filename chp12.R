############################################
###### chp12: Diagnosing Convergence #######
############################################

############################################
#### EX 12.1: Probit Model

# Data generation
genP <- function(n,beta,sigma,mu,tau){
  x <- rnorm(n,mean=mu,sd=tau)
  y <- runif(n)<pnorm(x*beta/sqrt(sigma))
  return(list(x=x,y=y))
}

# Posterior surface
likP <- function(x,y,lB,uB,lS,uS,prec=100){
  beta <- seq(lB,uB,(uB-lB)/(prec-1))[1:prec]
  sigma <- seq(lS,uS,(uS-lS)/(prec+9))[1:(prec+10)]
  div <- beta%*%t(1/sqrt(sigma))
  sur <- matrix(0,nrow=prec,ncol=(prec+10))
  for (i in 1:prec){
    for (j in 1:(prec+10)){
      a <- x*div[i,j]
      sur[i,j] <- sum( pnorm(a,log=T)*y + pnorm(-a,log=T)*(1-y) )
    }
  }
  
  # prior: exp(1) on 1/sigma and N(0,25) on beta
  for (j in 1:(prec+10)){
    sur[,j] <- sur[,j]-2*log(sigma[j])-1/sigma[j]}
  
  for (i in 1:prec){
    sur[i,] <- sur[i,]-beta[i]^2/50
  }
  return(list(l=sur,b=beta,s=sigma))
}

# posterior
lik1 <- function(x,y,B,S)
  sum(pnorm(x*B/sqrt(S),log=T)*y +
        pnorm(-x*B/sqrt(S),log=T)*(1-y) ) -
  2*log(S) - B^2/50 - 1/S

#Gibbs sampling
gibP <- function(m,x,y,B0,S0){
  beta <- 0*(1:m)
  sigma <- beta
  beta[1] <- B0
  sigma[1] <- S0
  like <- 1:m
  oldL <- lik1(x,y,B0,S0)
  like[1] <- oldL
  for (i in 2:m){
    
    # 1. Simulation of beta
    propB <- beta[i-1]+rnorm(1)
    propL <- lik1(x,y,propB,sigma[i-1])
    if (log(runif(1))<propL-oldL){
      beta[i] <- propB
      oldL <- propL
    }
    else{
      beta[i] <- beta[i-1]
    }
    # 2. Simulation of sigma
    propS <- exp(log(sigma[i-1])+.2*rnorm(1))
    propL <- lik1(x,y,beta[i],propS)
    if (log(runif(1))<propL-oldL-log(sigma[i-1])+log(propS)){
      sigma[i] <- propS
      oldL <- propL
    }
    else{
      sigma[i] <- sigma[i-1]
    }
    like[i] <- oldL
  }
  list(b=beta,s=sigma,l=like)
}

############################################
#### EX 12.7: Bimodal target
set.seed(2)

## Density of interest : exp(-x^2/2) (4*(x-.3)^2+.01)
## with normalizing constant 1 / (sqrt(2*pi)*(.01+4*(1+.3^2)))
f <- function(x){
  exp(-x^2/2) * (4*(x-.3)^2+.01) / (sqrt(2*pi)*(.01+4*(1+.3^2)))
}

## Generation of a slow/fast random walk
thet <- 0*(1:2000)
for (i in 2:2000){
  prop <- thet[i-1]+.3*rnorm(1)
  if (runif(1)<f(prop)/f(thet[i-1]))
    thet[i] <- prop
  else
    thet[i] <- thet[i-1]
}

## Corresponding sequence of masses
mass <- 0*(1:2000)
that <- thet[1]
for (i in 2:2000){
  that <- sort(c(that,thet[i]))
  mass[i] <- sum((that[2:i]-that[1:(i-1)])*f(that[1:(i-1)]))
}

## Plots
plot(mass,type="l",ylab="mass",col="sienna4",lwd=2)
par(new=T);plot(thet,pch=5,axes=F,cex=.3,col="steelblue2",ylab="")





