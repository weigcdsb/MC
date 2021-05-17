############################################
##### chp2: random variable generation #####
############################################


############################################
#### EX 2.11: Beta generation (Johnk's Thm)
set.seed(1)

n <- 10000
U <- runif(n)
V <- runif(n)

## generation
alph <- 2
bet <- 3
beta_samp <- 
ifelse((U^(1/alph) + V^(1/bet)) <= 1,
       (U^(1/alph))/(U^(1/alph) + V^(1/bet)), NA)

x <- seq(0, 1, length.out = 1000)
fx <- dbeta(x, alph, bet)

hist(beta_samp, freq = F)
lines(x, fx, col = 'red', lwd = 2)

## acceptance rate
alpha <- .1*(1:100)
probs <- rep(0,100)
for (i in 1:100){
  probs[i] <- length( U[(U^(1/alpha[i])+V^(1/alpha[i]))<=1] )/n
}

plot(alpha,probs,type="l",xlab=expression(alpha==beta),
     ylab="accept freq",lwd=2)

############################################
#### EX 2.16: Beta simulation (accept-reject)
#### Beta(2.7,6.3) based on uniforms
set.seed(2)
n <- 5000

ys <- runif(n)
maxd <- dbeta(1.7/(1.7+5.3),2.7,6.3)
us <- maxd*runif(n)
val <- seq(0,1.,.01)
valf <- dbeta(val,shape1=2.7,shape2=6.3)

as <- ifelse(us < dbeta(ys,2.7,6.3), us, NA)
plot(ys, as)
lines(val, valf, lwd = 2, col = 'red')


plot(val,valf,type="l",xlab="x",ylab="density",lwd=2)
polygon(c(val,rev(val)),c(valf,0*valf),col="cyan")
points(ys,us,cex=.4,pch=20)











