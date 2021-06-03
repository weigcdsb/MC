############################################
######### chp13: Perfect Sampling ##########
############################################

############################################
#### EX 13.11: truncated normal distribution
set.seed(1)

## mu=(-2.4,1.8) and Sigma=((1,-2.2),(-2.2,4.4))
mu <- c(-2.4,1.8)
sigma <- matrix(c(1,-1.2,-1.2,4.4),nrow=2,ncol=2)
sigmainv <- chol2inv(chol(sigma)) # choleski decomposition
# sigmainv <- solve(sigma)

## density
densi <- function(a,b){
  exp(-(c(a/(1-a),b/(1-b))-mu)%*%sigmainv%*%(c(a/(1-a),b/(1-b))-mu)/2)/((1-b)*(1-a))^2
}

## computation of the levels sets
x <- seq(0,3,.01)
y <- seq(0,10,.01)
zx <- x%*%t(rep(1,length(y)))
zy <- rep(1,length(x))%*%t(y)

## level plot
level <- sigmainv[1,1]*(zx-mu[1])^2+2*sigmainv[1,2]*(zx-mu[1])*(zy-mu[2])+sigmainv[2,2]*(zy-mu[2])^2
zx <- 1/(1+zx)
zy <- 1/(1+zy)
level <- exp(-0.5*level)/(zx*zy)^2
zx <- x/(1+x)
zy <- y/(1+y)
image(zx,zy,level,xlab=expression(x/1+x),ylab=expression(y/1+y))
contour(zx,zy,level,add=T)







