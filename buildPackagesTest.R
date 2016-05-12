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

library(mice)
library(statmod)
library(MCMCpack)
library(mitools)
library(mvtnorm)
library(mibrr)

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

nObs <- 125
nPreds <- 100
r2 <- 0.5
collin <- 0.3
betaRange <- c(0.2, 0.6)
meanRange <- c(0, 1)

dat1 <- simulateData(nObs = nObs,
                     nPreds = nPreds,
                     r2 = r2,
                     collin = collin,
                     betaRange = betaRange,
                     meanRange = meanRange)

testOut <- ben(rawData = dat1,
               y       = "y",
               X       = paste0("x", c(1 : nPreds))
               )

testOut2 <- bl(rawData = dat1,
               y       = "y",
               X       = paste0("x", c(1 : nPreds))
               )

lmOut <- lm(as.matrix(dat1[ , 1]) ~ as.matrix(dat1[ , -1]))

nTests <- 10
mseMat <- matrix(NA, nTests, 3)
for(i in 1 : nTests) {
    testDat <- simulateData(nObs = nObs,
                            nPreds = nPreds,
                            r2 = r2,
                            collin = collin,
                            betaRange = betaRange,
                            meanRange = meanRange)
    
    predOut1 <- predictMibrr(object  = testOut,
                             newData = as.matrix(testDat[ , -1]),
                             nDraws  = 500)
    
    predOut2 <- predictMibrr(object  = testOut2,
                             newData = as.matrix(testDat[ , -1]),
                             nDraws  = 500)
    
    predOut3 <- predict.lm(lmOut, newdata = testDat[ , -1])
    
    mseMat[i, 1] <- mean((rowMeans(predOut1) - dat1$y)^2)
    mseMat[i, 2] <- mean((rowMeans(predOut2) - dat1$y)^2)
    mseMat[i, 3] <- mean((predOut3 - dat1$y)^2)
}

colMeans(mseMat)

dens0 <- density(dat3$y)
dens1 <- density(rowMeans(predOut1))
dens2 <- density(rowMeans(predOut2))
dens3 <- density(predOut3)

yLim <- c(0, max(dens0$y, dens1$y, dens2$y, dens3$y))

plot(dens0, ylim = yLim)
lines(dens1, col = "red")
lines(dens2, col = "blue")
lines(dens3, col = "green")

