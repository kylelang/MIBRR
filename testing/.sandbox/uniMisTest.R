### Title:    Testing Univariate Missing Patterns Stuff
### Author:   Kyle M. Lang
### Created:  2018-NOV-19
### Modified: 2018-NOV-19

library(MIBRR)

nVars <- 65
nObs  <- 100

sigma       <- matrix(0.3, nVars, nVars)
diag(sigma) <- 1.0

X           <- rmvnorm(nObs, rep(1, nVars), sigma)
colnames(X) <- paste0("x", 1 : nVars)

rMat <- matrix(as.logical(rbinom(4 * nObs, 1, 0.3)), ncol = 4)

X[, 1 : 4][rMat] <- NA

out <- vanilla(data       = as.data.frame(X[ , 1 : 51]),
               targetVars = paste0("x", 1 : 4),
               ridge      = 0.75)

out$gibbsOut$x2$imps
