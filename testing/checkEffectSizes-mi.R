### Title:    Run Tests to Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2019-01-01
### Modified: 2019-11-25

rm(list = ls(all = TRUE))

expNum <- 1

nReps <- 500

sparse      <- TRUE
nImps       <- 100
nObs        <- 100
startRep    <- 1
stopRep     <- 16
clusterSize <- 8
outDir      <- NULL
verbose     <- FALSE
zCor        <- 0.5

resDir <- "output"

source("initScript-mcem.R")

nVars <- with(parms, length(c(y, X))) + control$nAux

r2    <- rep(NA, nReps)
coefs <- matrix(NA, nReps, nVars)
for(i in 1 : nReps) {
    dat <- as.data.frame(scale(simData(parms = parms, control = control)))
    out <- lm(y ~ ., data = dat)

    r2[i]      <- summary(out)$r.squared
    coefs[i, ] <- coef(out)
}

colMeans(coefs)
apply(coefs, 2, median)

mean(r2)
median(r2)

hist(r2)
median(r2)
