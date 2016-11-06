### Title:    Build R Packages
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2016-NOV-05
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))

#install.packages(c("statmod", "MCMCpack"),
#                 repos = "http://rweb.quant.ku.edu/cran")

library(RcppEigen)

system("rm source/mibrr/src/RcppExports.cpp \
        rm source/mibrr/R/RcppExports.R \
        rm source/mibrr/src/*.o source/mibrr/src/*.so")
Rcpp::compileAttributes("source/mibrr")
install.packages("source/mibrr", repos = NULL, type = "source")

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

miceOut <- mice(mibrrExampleData,
                m = 1,
                maxit = 100)

dat0 <- complete(miceOut)
dat1 <- data.frame(mibrrExampleData[ , 1 : 5], dat0[ , 6 : 17])

debug(miben)
undebug(miben)

testOut <- miben(rawData      = dat1,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                                        #iterations   = c(2, 1),
                                        #sampleSizes  = c(100, 500, 1000),
                 returnParams = TRUE,
                 verbose      = TRUE,
                 control      = list(simpleIntercept = TRUE,
                                     adaptScales     = TRUE)
                 )

fitOut <- lapply(testOut$imps,
                 FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                 )
MIcombine(fitOut)

summary(lm(y ~ x1 + x2 + x3, data = dat1))


testOut2 <- mibl(rawData      = mibrrExampleData,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                 returnParams = TRUE)

fitOut2 <- lapply(testOut2$imps,
                  FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                  )
MIcombine(fitOut2)


## Test BEN and BL:
library(mibrr)

alpha  <- 0.5
nPreds <- 175

parms             <- list()
parms$nObs        <- 100
parms$nPreds      <- nPreds
parms$r2          <- 0.2
parms$collin      <- 0.3
parms$beta        <- matrix(c(alpha, runif(nPreds, 0.3, 0.6)))
parms$means       <- runif(nPreds, 0, 1)
parms$scales      <- rep(1, nPreds)
parms$center      <- TRUE
parms$scale       <- TRUE
parms$simpleInt   <- FALSE
parms$verbose     <- FALSE
parms$postN       <- 5000
parms$adaptScales <- TRUE

testFun(1, parms)

rp <- 1

colMeans(dat1)

testFun <- function(rp, parms) {
    nObs   <- parms$nObs
    nPreds <- parms$nPreds
    r2     <- parms$r2
    collin <- parms$collin
    beta   <- parms$beta
    means  <- parms$means
    scales <- parms$scales

    dat1 <- simulateData(nObs   = nObs,
                         nPreds = nPreds,
                         r2     = r2,
                         collin = collin,
                         beta   = beta,
                         means  = means,
                         scales = scales)
    
    testOut <- ben(rawData = dat1,
                   y       = "y",
                   X       = paste0("x", c(1 : nPreds)),
                   verbose = parms$verbose,
                   control =
                       list(center          = parms$center,
                            scale           = parms$scale,
                            adaptScales     = parms$adaptScales,
                            simpleIntercept = parms$simpleInt)
                   )
    
    testOut2 <- bl(rawData = dat1,
                   y       = "y",
                   X       = paste0("x", c(1 : nPreds)),
                   verbose = parms$verbose,
                   control =
                       list(center          = parms$center,
                            scale           = parms$scale,
                            adaptScales     = parms$adaptScales,
                            simpleIntercept = parms$simpleInt)
                   )
    
    testForm <-
        as.formula(paste0("y ~ ", paste0("x", c(1 : nPreds), collapse = " + ")))
    lmOut <- lm(testForm, data = dat1)
    
    nTests <- 100
    mseMat <- matrix(NA, nTests, 3)
    for(i in 1 : nTests) {
        testDat <- simulateData(nObs   = nObs,
                                nPreds = nPreds,
                                r2     = r2,
                                collin = collin,
                                beta   = beta,
                                means  = means,
                                scales = scales)
        
        predOut1 <- predictMibrr(object  = testOut,
                                 newData = as.matrix(testDat[ , -1])
                                 )
        
        predOut2 <- predictMibrr(object  = testOut2,
                                 newData = as.matrix(testDat[ , -1])
                                 )
        
        predOut3 <- predict(lmOut, newdata = testDat[ , -1])
        
        mseMat[i, 1] <- mean((predOut1 - dat1$y)^2)
        mseMat[i, 2] <- mean((predOut2 - dat1$y)^2)
        mseMat[i, 3] <- mean((predOut3 - dat1$y)^2)
    }
    
    outMat           <- rbind(colMeans(mseMat), apply(mseMat, 2, var))
    colnames(outMat) <- c("MIBEN", "MIBL", "MLR")
    rownames(outMat) <- c("MSE", "var(MSE)")
    outMat
}# END testFun()

library(parallel)
nReps <- 10
mseList <- mclapply(c(1 : nReps),
                    FUN = testFun,
                    parms = parms,
                    mc.cores = 4)

mseList

?predict.lm


testFun(2, parms)
?predict.lm

dens0 <- density(dat3$y)
dens1 <- density(rowMeans(predOut1))
dens2 <- density(rowMeans(predOut2))
dens3 <- density(predOut3)

yLim <- c(0, max(dens0$y, dens1$y, dens2$y, dens3$y))

plot(dens0, ylim = yLim)
lines(dens1, col = "red")
lines(dens2, col = "blue")
lines(dens3, col = "green")

ls(testOut)

plot(testOut2$lambdaHistory[[1]][ , 1], type = "l")









mod1 <- paste(
    paste0("F",
           colnames(dat1),
           " =~ 1*",
           colnames(dat1),
           "\n"),
    collapse = "")

cat(mod1)

## Estimate the sufficient statistics with FIML:
out1 <- lavaan(model           = mod1,
               data            = dat1,
               int.ov.free     = FALSE,
               int.lv.free     = TRUE,
               auto.var        = TRUE,
               auto.fix.single = TRUE,
               missing         = "fiml")
        
## Store the item means and scales:
dataMeans  <- as.vector(inspect(out1, "coef")$alpha)
dataScales <- sqrt(diag(inspect(out1, "coef")$psi))

(dataMeans - colMeans(dat1)) / colMeans(dat1)
(dataScales - apply(dat1, 2, sd)) / apply(dat1, 2, sd)

names(env$dataMeans) <- names(env$dataScales) <- colnames(env$rawData)
