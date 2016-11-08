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

library(mibrr)
library(mitools)

## Test MIBEN and MIBL:
data(mibrrExampleData)

miceOut <- mice(mibrrExampleData,
                m = 1,
                maxit = 100)
dat1 <- complete(miceOut)
dat2 <- dat1
dat2$y <- mibrrExampleData$y

head(dat2)
debug(miben)
undebug(miben)

outList <- list()

testOut <- miben(data         = mibrrExampleData,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                 returnParams = TRUE,
                 verbose      = TRUE,
                 control      =
                     list(fimlStarts      = FALSE,
                          simpleIntercept = TRUE,
                          adaptScales     = TRUE)
                 )

colMeans(is.na(testOut$imps[[1]]))

colMeans(is.na(do.call(rbind, testOut$imps)))

testOut <- mibl(data         = mibrrExampleData,
                targetVars   = c("y", paste0("x", c(1 : 3))),
                ignoreVars   = "idNum",
                returnParams = TRUE,
                verbose      = TRUE,
                control      =
                    list(fimlStarts      = TRUE,
                         simpleIntercept = TRUE,
                         adaptScales     = TRUE)
                )

fitOut <- lapply(testOut$imps,
                 FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                 )
MIcombine(fitOut)


outList

round(coef(tmp), 3)

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
