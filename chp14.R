####################################################
# chp14: iterated & sequential importance sampling #
####################################################

############################################
#### EX 14.2: target tracking
set.seed(1)

datasim <- function(size,tau,eta){
  #Initial location & speed
  x <- 0:size
  x[1] <- 0.01
  xp <- 1.2
  y <- 0:size
  y[1] <- 5
  yp <- .3
  z <- 0:size
  z[1] <- (y[1]/x[1])+eta*rnorm(1)
  for (i in 2:(size+1)){
    xp <- xp+tau*rnorm(1)
    x[i] <- x[i-1]+xp
    yp <- yp+tau*rnorm(1)
    y[i] <- y[i-1]+yp
    z[i] <- (y[i]/x[i])+eta*rnorm(1)
  }
  list(angle=z,x=x,y=y)
}

dataplot <- function(x,y,z){
  # Plots the tracks of the target and the
  # corresponding angles
  plot(x,y,type="l",xlab="x",ylab="y",
       xlim=c(min(0,x),max(x)+1),ylim=c(min(0,y),max(y)+1))
  points(x,y,pch="+",cex=1.5)
  for (i in 1:101)
    lines(c(0,10*max(x,y)),c(0,10*max(x,y)*z[i]),lty=3)
}

sis <- function(size,npart,z,tau,eta){
  # Simulate a particle population for size
  # iterations and npart particles w/o resampling
  # Initialisation
  omega <- matrix(0,ncol=npart,nrow=size)
  xpart <- omega
  ypart <- omega
  # Run
  vx <- 1.2+tau*rnorm(npart)
  vy <- .3+tau*rnorm(npart)
  xpart[1,] <- 1.2+vx
  ypart[1,] <-  5 +vy
  omega[1,] <- (z[2]-ypart[1,]/xpart[1,])^2
  for (t in 2:size){
    vx <- vx+tau*rnorm(npart)
    vy <- vy+tau*rnorm(npart)
    xpart[t,] <- xpart[t-1,]+vx
    ypart[t,] <- ypart[t-1,]+vy
    omega[t,] <- (z[t+1]-ypart[t,]/xpart[t,])^2
  }
  omega <-- -apply(omega,2,cumsum)/(2.*eta^2) #cumulated sums
  best <- (t(apply(t(omega),2,order)))[,npart]  # should give the best at each iteration
  bestun <- best+npart*((1:size)-1)
  bestend <- best[size]
  list(weights=t(omega)[bestun],x=t(xpart)[bestun],y=t(ypart)[bestun],
       xend=xpart[,bestend],yend=ypart[,bestend])
}

ris <- function(size,npart,z,tau,eta){
  # Simulate a particle population for size
  # iterations and npart particles w resampling
  # Initialisation
  omega <- matrix(0,ncol=npart,nrow=size)
  xpart <- omega
  xstep <- xpart
  ypart <- omega
  ystep <- ypart
  # Run : starts from true location and speed
  vx <- 1.2+tau*rnorm(npart)
  vy <- .3+tau*rnorm(npart)
  xpart[1,] <- .01+vx
  ypart[1,] <- 5+vy
  omega[1,] <- exp(-(z[2]-ypart[1,]/xpart[1,])^2/(2*eta^2))
  omega[1,] <- omega[1,]/sum(omega[1,])
  survl <- sample(1:npart,npart,replace=TRUE,prob=omega[1,])
  x <- xpart[1,survl]
  y <- ypart[1,survl]
  xpart[1,] <- x
  ypart[1,] <- y
  xstep[1,] <- x
  ystep[1,] <- y
  for (t in 2:size){
    vx <- vx+tau*rnorm(npart)
    vy <- vy+tau*rnorm(npart)
    xpart[t,] <- xpart[t-1,]+vx
    ypart[t,] <- ypart[t-1,]+vy
    omega[t,] <- exp(-(z[t+1]-ypart[t,]/xpart[t,])^2/(2*eta^2))
    omega[t,] <- omega[t,]/sum(omega[t,])
    survl <- sample(1:npart,npart,replace=TRUE,prob=omega[t,])
    for (j in 1:t){ # Reallocate the whole path
      x <- xpart[j,survl]
      y <- ypart[j,survl]
      xpart[j,] <- x
      ypart[j,] <- y
    }
    xstep[t,] <- xpart[t,] # Only change the current particles
    ystep[t,] <- ypart[t,]
  }
  list(weights=omega,x=xpart,y=ypart,xst=xstep,yst=ystep)
}

whole <- function(size,npart,tau,eta){
  # Create the data
  res <- datasim(size,tau,eta)
  # Runs the SIS procedure
  ras <- sis(size,npart,res$angle,tau,eta)
  # Plots the results
  plot(res$x,res$y,type="l",xlab="x",ylab="y",lwd=2,
       xlim=c(min(0,res$x),max(res$x)+1),ylim=c(min(0,res$y),max(res$y)+1))
  lines(ras$x,ras$y,lty=2,col="sienna4",lwd=2)
  lines(ras$xend,ras$yend,lty=3,col="steelblue",lwd=2)
  title(main="SIS")
  # Runs the RIS procedure
  rus <- ris(size,npart,res$angle,tau,eta)
  # Plots the results
  plot(res$x,res$y,type="l",xlab="x",ylab="y",lwd=2,
       xlim=c(min(0,res$x),max(res$x)+1),ylim=c(min(0,res$y),max(res$y)+1))
  lines(apply(rus$x,1,mean),apply(rus$y,1,mean),lty=2,col="sienna4",lwd=2)
  lines(apply(rus$xst,1,mean),apply(rus$yst,1,mean),lty=3,col="steelblue",lwd=2)
  title(main="RIS")
  # Ranges
  plot(res$x,res$y,type="l",xlab="x",ylab="y",lwd=2,
       xlim=c(min(0,res$x),max(res$x)+1),ylim=c(min(0,res$y),max(res$y)+1))
  sampin <- sample(1:npart,100)
  for (i in sampin)
    lines(rus$x[,i],rus$y[,i],col="gray")
  lines(res$x,res$y,lwd=2)
  title(main="Final RIS range")
  plot(res$x,res$y,type="l",xlab="x",ylab="y",lwd=2,
       xlim=c(min(0,res$x),max(res$x)+1),ylim=c(min(0,res$y),max(res$y)+1))
  sampin_sample(1:npart,100)
  for (i in sampin)
    lines(rus$xst[,i],rus$yst[,i],col="gray")
  lines(res$x,res$y,lwd=2)
  title(main="Stepwise RIS range")
}