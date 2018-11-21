### Title:    Testing Univariate Missing Patterns Stuff
### Author:   Kyle M. Lang
### Created:  2018-NOV-19
### Modified: 2018-NOV-19

library(MIBRR)

sigma       <- matrix(0.3, 5, 5)
diag(sigma) <- 1.0

X           <- rmvnorm(1000, rep(1, 5), sigma)
colnames(X) <- paste0("x", 1 : 5)

X[as.logical(rbinom(1000, 1, 0.3)), 2] <- NA

out <- miben(data = as.data.frame(X), targetVars = "x2")

out$gibbsOut$x2$imps
