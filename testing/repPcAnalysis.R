### Title:    Run Tests to Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2019-01-01
### Modified: 2019-11-28

rm(list = ls(all = TRUE))
    
library(parallel)
library(lars)
library(MIBRR)

data(diabetes)

y <- diabetes$y
X <- diabetes$x

dat1 <- data.frame(cbind(y, X))

fit      <- lm(y ~ X)
lamStart <- ncol(X) * summary(fit)$sigma / sum(abs(coef(fit)[-1]))

out <- list()
for(rp in 1 : 4) 
    out[[rp]] <- bl(data        = dat1,
                    y           = "y",
                    X           = colnames(X),
                    iterations  = c(500, 100),
                    sampleSizes = list(rep(250, 2), rep(500, 2), rep(1000, 2)),
                    control     = list(
                        lambda1Starts =
                            lamStart + abs(rnorm(1, 0, sqrt(lamStart)))
                    )
                    )

getParams(out, "y")

out$y$lambdaHistory

plot(out[[1]]$lambdaHistory$y[ , 1], type = "l")
lines(out[[2]]$lambdaHistory$y[ , 1], col = "red")
lines(out[[3]]$lambdaHistory$y[ , 1], col = "blue")
lines(out[[4]]$lambdaHistory$y[ , 1], col = "green")

dat2 <- as.data.frame(scale(dat1))

fit      <- lm(y ~ ., data = dat2)
lamStart <- ncol(X) * summary(fit)$sigma / sum(abs(coef(fit)[-1]))

out2 <- list()
for(rp in 1 : 4) 
    out2[[rp]] <- bl(data        = dat2,
                     y           = "y",
                     X           = colnames(X),
                     iterations  = c(500, 100),
                     sampleSizes = list(rep(250, 2), rep(500, 2), rep(1000, 2)),
                     control     = list(
                         lambda1Starts =
                             lamStart + abs(rnorm(1, 0, sqrt(lamStart)))
                     )
                     )

plot(out2[[1]]$lambdaHistory$y[ , 1], type = "l")
lines(out2[[2]]$lambdaHistory$y[ , 1], col = "red")
lines(out2[[3]]$lambdaHistory$y[ , 1], col = "blue")
lines(out2[[4]]$lambdaHistory$y[ , 1], col = "green")

pars1 <- getParams(out[[1]], "y")
pars2 <- getParams(out[[2]], "y")

pars1$beta
pars2$beta

med1 <- apply(pars1$beta, 2, median)
med2 <- apply(pars2$beta, 2, median)

lmL1 <- sum(abs(coef(fit)[-1]))
blL1 <- sum(abs(med2[-1]))

lmL1
blL1

+blL1 / lmL1

plot(med1)
med1
coef(fit)
