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

## Test MIBEN and MIBL:
data(mibrrExampleData)

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

testFun(1, parms)

rp <- 1


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

    rMat <- matrix(as.logical(rbinom(prod(dim(dat1)), 1, .2)), nrow(dat1))
    dat1[rMat] <- NA
    
    testOut <- ben(data    = dat1,
                   y       = "y",
                   X       = paste0("x", c(1 : nPreds)),
                   verbose = parms$verbose,
                   control =
                       list(adaptScales     = parms$adaptScales,
                            simpleIntercept = parms$simpleInt)
                   )
    
    testOut2 <- bl(data    = dat1,
                   y       = "y",
                   X       = paste0("x", c(1 : nPreds)),
                   verbose = parms$verbose,
                   control =
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
