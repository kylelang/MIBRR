
x <- matrix(runif(100))

crossprod(x - mean(x)) / 99

var(x)

library(optimx)

?optimx

## Show multiple outputs of optimx using all.methods
                                        # genrose function code
genrose.f<- function(x, gs=NULL){ # objective function
    ## One generalization of the Rosenbrock banana valley function (n parameters)
    n <- length(x)
    if(is.null(gs)) { gs=100.0 }
    fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
    return(fval)
}

genrose.g <- function(x, gs=NULL){
                                        # vectorized gradient for genrose.f
                                        # Ravi Varadhan 2009-04-03
    n <- length(x)
    if(is.null(gs)) { gs=100.0 }
    gg <- as.vector(rep(0, n))
    tn <- 2:n
    tn1 <- tn - 1
    z1 <- x[tn] - x[tn1]^2
    z2 <- 1 - x[tn]
    gg[tn] <- 2 * (gs * z1 - z2)
    gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
    return(gg)
}

genrose.h <- function(x, gs=NULL) { ## compute Hessian
    if(is.null(gs)) { gs=100.0 }
    n <- length(x)
    hh<-matrix(rep(0, n*n),n,n)
    for (i in 2:n) {
        z1<-x[i]-x[i-1]*x[i-1]
        z2<-1.0-x[i]
        hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
        hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
        hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
        hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
    }
    return(hh)
}

startx<-4*seq(1:10)/3.
ans8<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, 
             control=list(all.methods=TRUE, save.failures=TRUE, trace=0, kkt = FALSE), gs=10)

ans8

coef(ans8)
ans8[ , c("convcode", "kkt1", "kkt2")]


x <- rnorm(1000, 5, 25)

y <- 2 + 3 * x

beta1 <- matrix(coef(lm(y ~ x)))
yHat1 <- cbind(1, x) %*% beta1

beta2 <- matrix(coef(lm(y ~ scale(x))))
yHat2 <- cbind(1, scale(x)) %*% beta2

sum(yHat1 - yHat2)
