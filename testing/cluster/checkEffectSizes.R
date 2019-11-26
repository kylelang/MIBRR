### Title:    Run Tests to Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2019-01-01
### Modified: 2019-11-25

rm(list = ls(all = TRUE))
    
library(parallel)

verbose <- FALSE
resDir  <- "output/"
nCores  <- 2
nObs    <- 100
nVars   <- 10
pm      <- 0.0
xCor    <- 0.0
nReps   <- 500
mi      <- FALSE

source("initScript-simple.R")

r2    <- rep(NA, nReps)
coefs <- matrix(NA, nReps, nVars)
for(i in 1 : nReps) {
    dat <- as.data.frame(scale(genSimpleData(parms = parms, pm = pm)))
    out <- lm(x1 ~ ., data = dat)

    r2[i]      <- summary(out)$r.squared
    coefs[i, ] <- coef(out)
}

colMeans(coefs)
apply(coefs, 2, median)

mean(r2)
median(r2)

hist(r2)
median(r2)
