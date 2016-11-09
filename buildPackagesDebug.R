### Title:    Build R Packages
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2016-NOV-05
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))

library(RcppEigen)

system("rm source/mibrr/src/RcppExports.cpp \
        rm source/mibrr/R/RcppExports.R \
        rm source/mibrr/src/*.o source/mibrr/src/*.so")
Rcpp::compileAttributes("source/mibrr")
install.packages("source/mibrr", repos = NULL, type = "source")

library(mibrr)
library(mitools)

data(mibrrExampleData)

## Use mice() as a baseline comparison:
miceOut <- mice(mibrrExampleData,
                m = 100,
                maxit = 20)
miceImps <- list()
for(m in 1 : 100) miceImps[[m]] <- complete(miceOut, m)

miceFits <- lapply(miceImps,
                   FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                   )
micePooled <- MIcombine(miceFits)
micePooled

## Test MIBEN and MIBL:
mibenOut <- miben(data           = mibrrExampleData,
                  targetVars     = c("y", paste0("x", c(1 : 3))),
                  ignoreVars     = "idNum",
                  returnConvInfo = TRUE,
                  returnParams   = TRUE,
                  verbose        = TRUE,
                  control        =
                      list(fimlStarts      = TRUE,
                           simpleIntercept = TRUE,
                           adaptScales     = TRUE)
                  )

mibenFits <- lapply(mibenOut$imps,
                    FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                    )
mibenPooled <- MIcombine(mibenFits)
mibenPooled

par(mfcol = c(2, 4))
for(v in 1 : 4) {
    tmp <- mibenOut$lambdaHistory[[v]]

    plot(tmp[ , 1], type = "l", main = "Lambda1")
    plot(tmp[ , 2], type = "l", main = "Lambda2")
}

any(unlist(mibenOut$rHats) > 1.1)
colMeans(is.na(do.call(rbind, mibenOut$imps)))

miblOut <- mibl(data           = mibrrExampleData,
                targetVars     = c("y", paste0("x", c(1 : 3))),
                ignoreVars     = "idNum",
                returnConvInfo = TRUE,
                returnParams   = TRUE,
                verbose        = TRUE,
                control        =
                    list(fimlStarts      = TRUE,
                         simpleIntercept = TRUE,
                         adaptScales     = TRUE)
                )

miblFits <- lapply(miblOut$imps,
                   FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                   )
miblPooled <- MIcombine(miblFits)
miblPooled

par(mfcol = c(2, 2))
for(v in 1 : 4) {
    tmp <- miblOut$lambdaHistory[[v]]
    plot(tmp, type = "l", main = "Lambda")
}

any(unlist(mibenOut$rHats) > 1.1)
colMeans(is.na(do.call(rbind, mibenOut$imps)))

## Test BEN and BL:

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
parms$verbose     <- TRUE
parms$postN       <- 5000
parms$adaptScales <- TRUE
parms$iterations  <- c(2000, 25)

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
      
    testOut <- ben(data       = dat1,
                   y          = "y",
                   X          = paste0("x", c(1 : nPreds)),
                   iterations = parms$iterations,
                   verbose    = parms$verbose,
                   control    =
                       list(adaptScales     = parms$adaptScales,
                            simpleIntercept = parms$simpleInt)
                   )
    
    testOut2 <- bl(data       = dat1,
                   y          = "y",
                   X          = paste0("x", c(1 : nPreds)),
                   iterations = parms$iterations,
                   verbose    = parms$verbose,
                   control    =
                       list(adaptScales     = parms$adaptScales,
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
        
        predOut1 <- predictMibrr(object    = testOut,
                                 newData   = as.matrix(testDat[ , -1]),
                                 targetVar = "y",
                                 nDraws    = 250)

        meanPred1 <- rowMeans(predOut1)
        varPred1  <- apply(predOut1, 1, var)
        
        predOut2 <- predictMibrr(object    = testOut2,
                                 newData   = as.matrix(testDat[ , -1]),
                                 targetVar = "y",
                                 nDraws    = 250)

        meanPred2 <- rowMeans(predOut2)
        varPred2  <- apply(predOut2, 1, var)
        
        predOut3 <- predict(lmOut, newdata = testDat[ , -1])
        
        mseMat[i, 1] <- mean((meanPred1 - dat1$y)^2)
        mseMat[i, 2] <- mean((meanPred2 - dat1$y)^2)
        mseMat[i, 3] <- mean((predOut3 - dat1$y)^2)
    }
    
    outMat           <- rbind(colMeans(mseMat), apply(mseMat, 2, var))
    colnames(outMat) <- c("MIBEN", "MIBL", "MLR")
    rownames(outMat) <- c("MSE", "var(MSE)")
    outMat
}# END testFun()

mean(varPred1)
mean(varPred2)

library(parallel)
nReps <- 10
mseList <- mclapply(c(1 : nReps),
                    FUN = testFun,
                    parms = parms,
                    mc.cores = 4)

tmp1 <- testOut$lambdaHistory$y
tmp2 <- testOut2$lambdaHistory$y

par(mfrow = c(1, 3))
plot(tmp1[ , "lambda1"], type = "l")
plot(tmp1[ , "lambda2"], type = "l")
plot(tmp2, type = "l")
