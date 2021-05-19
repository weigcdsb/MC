############################################
## chp4: controling Monte Carlo variance ###
############################################

############################################
#### EX 4.2: Cauchy prior (Cauchy(0, 1))
#### univariate monitoring
set.seed(1)

h1 <- function(x){x/(1+x^2)}
h2 <- function(x){1/(1+x^2)}

## Implementation of a real confidence band:
N <- 1000 # Total number of simulations
M <- 1000  # Number of parallel estimators
uper <- round(.95*M)
loer <- round(.05*M)
uppa <- 0*(1:N)
lowa <- uppa
avery <- lowa
var <- avery
var1 <- avery
var2 <- avery
cova <- avery
Xlo <- lowa
Xhi <- lowa
kapa <- lowa # normalizing constant
estim1 <- 0*(1:M)
estim2 <- 0*(1:M)
estim <- estim1
for (i in 1:N){
  norma <- rnorm(M)+2.5 # x 
  ash1 <- h1(norma)
  ash2 <- h2(norma)
  estim1 <- estim1+ash1
  estim2 <- estim2+ash2
  
  ## Empirical evaluation
  estim <- sort(estim1/estim2)
  avery[i] <- mean(estim)
  lowa[i] <- estim[loer]
  uppa[i] <- estim[uper]
  Xlo[i] <- min(estim)
  Xhi[i] <- max(estim)
  
  ## Asymptotic variance evaluation
  kapa[i] <- mean(estim2)/i
  var1[i] <- var(estim1)
  var2[i] <- var(estim2)
  cova[i] <- cov(estim1,estim2)
  var[i] <- (var1[i]-2*avery[i]*cova[i]+avery[i]^2*var2[i])/(i^2*kapa[i]^2)
}


## plot
plot(avery,ylim=c(.5,4.5),xlab="Iterations",ylab="",type="l",
     lwd=2.5,col="blue")
polygon(c(1:N,N:1),c(Xlo,Xhi[N:1]),col="grey")
lines(avery,lwd=2.5,col="blue")
lines(uppa,lty=2,lwd=2.5,col="blue")
lines(lowa,lty=2,lwd=2.5,col="blue")
lines(avery+sqrt(var),lty=5,lwd=2.5,col="red")
lines(avery-sqrt(var),lty=5,lwd=2.5,col="red")
avery[N]

############################################
#### EX 4.4: continuation of Ex 4.2
#### importance sampling version
set.seed(2)

h1 <- function(x){x*exp(-(x-2.5)^2/2)}
h2 <- function(x){exp(-(x-2.5)^2/2)}

## Implementation of a real confidence band:
N <- 1000 # Total number of simulations
M <- 5000  # Number of parallel estimators
uper <- round(.95*M)
loer <- round(.05*M)
uppa <- 0*(1:N)
lowa <- uppa
avery <- lowa
var <- avery
var1 <- avery
var2 <- avery
cova <- avery
Xlo <- lowa
Xhi <- lowa
kapa <- lowa # normalizing constant
estim1 <- 0*(1:M)
estim2 <- 0*(1:M)
estim <- estim1
for (i in 1:N){
  norma <- rnorm(M)/rnorm(M) # this is a Cauchy (0,1)
  ash1 <- h1(norma)
  ash2 <- h2(norma)
  estim1 <- estim1+ash1
  estim2 <- estim2+ash2
  
  ## Empirical evaluation
  estim <- sort(estim1/estim2)
  avery[i] <- mean(estim)
  lowa[i] <- estim[loer]
  uppa[i] <- estim[uper]
  Xlo[i] <- min(estim)
  Xhi[i] <- max(estim)
  
  ## Asymptotic variance evaluation
  kapa[i] <- mean(estim2)/i
  var1[i] <- var(estim1)
  var2[i] <- var(estim2)
  cova[i] <- cov(estim1,estim2)
  var[i] <- (var1[i]-2*avery[i]*cova[i]+avery[i]^2*var2[i])/(i^2*kapa[i]^2)
}

## plot
plot(avery,ylim=c(.5,4.5),xlab="Iterations",ylab="",
     type="l",lwd=2.5,col="blue")
polygon(c(1:N,N:1),c(Xlo,Xhi[N:1]),col="grey")
lines(avery,lwd=2.5,col="blue")
lines(uppa,lty=2,lwd=2.5,col="blue")
lines(lowa,lty=2,lwd=2.5,col="blue")
lines(avery+sqrt(var),lty=5,lwd=2.5,col="red")
lines(avery-sqrt(var),lty=5,lwd=2.5,col="red")
print(avery[N])












