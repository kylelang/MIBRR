### Title:    Build R Packages
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2016-MAY-05
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))

install.packages("statmod",
                 repos = "http://rweb.quant.ku.edu/cran")

library(RcppEigen)

system("rm -r packageSource/mibrr/src/nlopt/* \
        rm packageSource/mibrr/src/RcppExports.cpp \
        rm packageSource/mibrr/R/RcppExports.R \
        rm packageSource/mibrr/src/*.o packageSource/mibrr/src/*.so")
system("cp -r ~/data/software/miscPackages/nlopt-2.4.2/* packageSource/mibrr/src/nlopt/")
Rcpp::compileAttributes("packageSource/mibrr")
install.packages("packageSource/mibrr", repos = NULL, type = "source")

library(statmod)
library(MCMCpack)
library(mitools)
library(mvtnorm)
library(mibrr)
library(mice)

## Test MVN sampler:
nObs <- 500000
mvnMu <- rep(10, 3)
mvnSigma <- matrix(5, 3, 3)
diag(mvnSigma) <- 20

out1.1 <- mibrr::drawMVN(nObs, mvnMu, mvnSigma)
out1.2 <- rmvnorm(nObs, mvnMu, mvnSigma)

par(mfrow = c(1, 3))
for(i in 1 : length(mvnMu)) {
    plot(density(out1.1[ , i]), col = "red")
    lines(density(out1.2[ , i]), col = "blue")
}

## Test inverse gamma sampler:
gamShape <- 10
gamScale <- 10

out2.1 <- mibrr::drawInvGamma(nObs, gamShape, gamScale)
out2.2 <- rinvgamma(nObs, gamShape, gamScale)

plot(density(out2.1), col = "red")
lines(density(out2.2), col = "blue")

## Test inverse gaussian sampler:
igMu = 1
igLam = 2

out3.1 <- mibrr::drawInvGauss(nObs, igMu, igLam)
out3.2 <- rinvgauss(nObs, igMu, igLam)

plot(density(out3.1), col = "red")
lines(density(out3.2), col = "blue")

## Test the incomplte gamma calculation:
incGamShape <- 10
incGamCut <- 5

out4.1 <- mibrr::calcIncGamma(incGamShape, incGamCut, FALSE)
out4.2 <- pgamma(q = incGamCut,
                 shape = incGamShape,
                 lower = FALSE) * gamma(incGamShape)

out4.1 - out4.2

## Test MIBEN and MIBL:
data(mibrrExampleData)

debug(miben)
undebug(miben)

testOut <- miben(rawData      = mibrrExampleData,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                 returnParams = TRUE)

fitOut <- lapply(testOut$imps,
                 FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                 )
MIcombine(fitOut)

testOut2 <- mibl(rawData      = mibrrExampleData,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                 returnParams = TRUE)

fitOut2 <- lapply(testOut$imps,
                 FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                 )
MIcombine(fitOut2)


## Test BEN and BL:

dat1 <- complete(mice(mibrrExampleData, m = 1, maxit = 100), 1)

testOut <- ben(rawData      = dat1,
               y            = "y",
               X            = paste0("x", c(1 : 3)),
               returnParams = TRUE)

fitOut <- lapply(testOut$imps,
                 FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                 )
MIcombine(fitOut)

testOut2 <- mibl(rawData      = mibrrExampleData,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                 returnParams = TRUE)

fitOut2 <- lapply(testOut$imps,
                 FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                 )
MIcombine(fitOut2)
