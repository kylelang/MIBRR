### Title:    Build R Packages
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2016-NOV-05
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))

install.packages("enet", repos = "http://cloud.r-project.org")

library(RcppEigen)

system("rm source/mibrr/src/RcppExports.cpp \
        rm source/mibrr/R/RcppExports.R \
        rm source/mibrr/src/*.o source/mibrr/src/*.so")
Rcpp::compileAttributes("source/mibrr")
install.packages("source/mibrr", repos = NULL, type = "source")

library(mibrr)
library(mitools)
library(glmnet)

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

alpha          <- 50
nPreds         <- 10

parms                 <- list()
parms$nObs            <- 1000
parms$nPreds          <- nPreds
parms$r2              <- 0.5
parms$collin          <- 0.3
parms$beta            <- matrix(
    c(alpha, rep(0.3, ceiling(nPreds / 2)), rep(0, floor(nPreds / 2)))
)
parms$means           <- runif(nPreds, -1, 1)
parms$scales          <- rep(1, nPreds)
parms$simpleInt       <- TRUE
parms$verbose         <- FALSE
parms$adaptScales     <- TRUE
parms$iterations      <- c(2500, 25)
parms$latentStructure <- FALSE
parms$itemsPerFactor  <- 4
parms$itemReliability <- 0.8
parms$twoPhaseOpt     <- FALSE

testFun <- function(rp, parms) {
    nObs            <- parms$nObs
    nPreds          <- parms$nPreds
    r2              <- parms$r2
    collin          <- parms$collin
    beta            <- parms$beta
    means           <- parms$means
    scales          <- parms$scales
    latentStructure <- parms$latentStructure
    itemsPerFactor  <- parms$itemsPerFactor
    itemReliability <- parms$itemReliability
    
    dat1 <- simulateData(nObs            = nObs,
                         nPreds          = nPreds,
                         r2              = r2,
                         collin          = collin,
                         beta            = beta,
                         means           = means,
                         scales          = scales,
                         latentStructure = latentStructure,
                         itemsPerFactor  = itemsPerFactor,
                         itemReliability = itemReliability)
    
    nIVs <- ifelse(latentStructure, nPreds * itemsPerFactor, nPreds)

    X <- as.matrix(dat1[ , paste0("x", c(1 : nIVs))])
    y <- dat1$y
    
    benOut <- ben(data       = dat1,
                  y          = "y",
                  X          = paste0("x", c(1 : nIVs)),
                  iterations = parms$iterations,
                  verbose    = parms$verbose,
                  control    =
                      list(adaptScales     = parms$adaptScales,
                           simpleIntercept = parms$simpleInt,
                           twoPhaseOpt     = parms$twoPhaseOpt)
                  )
    
    blOut <- bl(data       = dat1,
                y          = "y",
                X          = paste0("x", c(1 : nIVs)),
                iterations = parms$iterations,
                verbose    = parms$verbose,
                control    =
                    list(adaptScales     = parms$adaptScales,
                         simpleIntercept = parms$simpleInt)
                )
    
    aVec <- seq(0, 1, 0.01)
    cvmVec <- lamVec <- rep(NA, length(aVec))
    for(i in 1 : length(aVec)) {
        cvOut <- cv.glmnet(x = X,
                           y = y,
                           nfolds = length(y),
                           alpha = aVec[i],
                           grouped = FALSE)
        
        best <- which(cvOut$lambda == cvOut$lambda.min)
        
        cvmVec[i] <- cvOut$cvm[best]
        lamVec[i] <- cvOut$lambda.min
    }
    
    alpha    <- aVec[which.min(cvmVec)]
    
    enOut    <- cv.glmnet(x       = X,
                          y       = y,
                          nfolds  = length(y),
                          alpha   = alpha,
                          grouped = FALSE)

    lassoOut <- cv.glmnet(x       = X,
                          y       = y,
                          nfolds  = length(y),
                          alpha   = 1,
                          grouped = FALSE)
    
    testForm <-
        as.formula(paste0("y ~ ", paste0("x", c(1 : nIVs), collapse = " + ")))

    lmOut <- lm(testForm, data = dat1)

    testDat <- simulateData(nObs            = nObs,
                            nPreds          = nPreds,
                            r2              = r2,
                            collin          = collin,
                            beta            = beta,
                            means           = means,
                            scales          = scales,
                            latentStructure = latentStructure,
                            itemsPerFactor  = itemsPerFactor,
                            itemReliability = itemReliability)
    
    newX <- as.matrix(testDat[ , paste0("x", c(1 : nIVs))])
    
    benPPD <- predictMibrr(object    = benOut,
                           newData   = as.matrix(testDat[ , -1]),
                           targetVar = "y",
                           nDraws    = 250)
    benPred <- rowMeans(benPPD)
    
    blPPD <- predictMibrr(object    = blOut,
                          newData   = as.matrix(testDat[ , -1]),
                          targetVar = "y",
                          nDraws    = 250)
    blPred <- rowMeans(blPPD)
    
    lmPred    <- predict(lmOut, newdata = testDat[ , -1])
    enPred    <- predict(enOut, newX)
    lassoPred <- predict(lassoOut, newX)
    
    outVec <- c(mean((benPred   - dat1$y)^2),
                mean((blPred    - dat1$y)^2),
                mean((enPred    - dat1$y)^2),
                mean((lassoPred - dat1$y)^2),
                mean((lmPred    - dat1$y)^2))
    names(outVec) <- c("MIBEN", "MIBL", "ENET", "LASSO", "MLR")
    outVec
}# END testFun()




testFun2 <- function(rp, parms) {
    nObs            <- parms$nObs
    nPreds          <- parms$nPreds
    r2              <- parms$r2
    collin          <- parms$collin
    beta            <- parms$beta
    means           <- parms$means
    scales          <- parms$scales
    latentStructure <- parms$latentStructure
    itemsPerFactor  <- parms$itemsPerFactor
    itemReliability <- parms$itemReliability
    
    dat1 <- simulateData(nObs            = nObs,
                         nPreds          = nPreds,
                         r2              = r2,
                         collin          = collin,
                         beta            = beta,
                         means           = means,
                         scales          = scales,
                         latentStructure = latentStructure,
                         itemsPerFactor  = itemsPerFactor,
                         itemReliability = itemReliability)
    
    nIVs <- ifelse(latentStructure, nPreds * itemsPerFactor, nPreds)

    outList <- list()
    outList$ben <- ben(data       = dat1,
                       y          = "y",
                       X          = paste0("x", c(1 : nIVs)),
                       iterations = parms$iterations,
                       verbose    = parms$verbose,
                       control    =
                           list(adaptScales     = parms$adaptScales,
                                simpleIntercept = parms$simpleInt,
                                twoPhaseOpt     = parms$twoPhaseOpt)
                       )
    
    outList$bl <- bl(data       = dat1,
                     y          = "y",
                     X          = paste0("x", c(1 : nIVs)),
                     iterations = parms$iterations,
                     verbose    = parms$verbose,
                     control    =
                         list(adaptScales     = parms$adaptScales,
                              simpleIntercept = parms$simpleInt)
                     )
    
    testForm <-
        as.formula(paste0("y ~ ", paste0("x", c(1 : nIVs), collapse = " + ")))

    tmpDat <- as.data.frame(scale(dat1, scale = FALSE))
    outList$lm <- lm(testForm, data = tmpDat)
    outList
}

?lm

out1 <- testFun2(1, parms)

benOut <- out1$ben
blOut  <- out1$bl
lmOut  <- out1$lm

benPPD <- predictMibrr(object    = benOut,
                       newData   = as.matrix(dat1[ , -1]),
                       targetVar = "y",
                       nDraws    = nrow(benOut$params$y$beta)
                       )

blPPD <- predictMibrr(object    = blOut,
                      newData   = as.matrix(dat1[ , -1]),
                      targetVar = "y",
                      nDraws    = nrow(benOut$params$y$beta)
                      )

lmPPD <- predict(lmOut, newData = dat1[ , -1])

benDens <- density(rowMeans(benPPD))
blDens  <- density(rowMeans(blPPD))
lmDens  <- density(lmPPD)
yDens   <- density(scale(dat1$y, scale = FALSE))

yLim <- range(benDens$y, blDens$y, lmDens$y, yDens$y)
xLim <- range(benDens$x, blDens$x, lmDens$x, yDens$x)

par(mfrow = c(1, 1))
plot(benDens, col = "red", ylim = yLim, xlim = xLim)
lines(blDens, col = "blue")
lines(lmDens, col = "green")
lines(yDens)

plot(enOut)

?plot.glmnet

?cv.glmnet
?glmnet
?enet

x2 <- apply(dat1, 2, as.numeric)

plot(lassoOut)
?plot.glmnet

lassoOut <- cv.glmnet(x = as.matrix(dat1[ , paste0("x", c(1 : nIVs))]),
                      y = dat1$y,
                      foldid = c(1 : nrow(dat1)),
                      alpha = 1)

ridgeOut <- cv.glmnet(x = as.matrix(dat1[ , paste0("x", c(1 : nIVs))]),
                      y = dat1$y,
                      foldid = c(1 : nrow(dat1)),
                      alpha = 0)

X <- as.matrix(dat1[ , paste0("x", c(1 : nIVs))])
y <- dat1$y


(is.na(cvmVec))
lamVec

warnings()

ls(ridgeOut)

ridgeOut$glmnet.fit
ridgeOut$lambda.min
lassoOut$lambda.min

library(enet)

out1 <- testFun(1, parms)
out1

mean(varPred1)
mean(varPred2)

library(parallel)
nReps <- 10
mseList <- mclapply(c(1 : nReps),
                    FUN = testFun,
                    parms = parms,
                    mc.cores = 4)

tmp1 <- benOut$lambdaHistory$y
tmp2 <- blOut$lambdaHistory$y

par(mfrow = c(1, 3))
plot(tmp1[ , "lambda1"], type = "l")
plot(tmp1[ , "lambda2"], type = "l")
plot(tmp2, type = "l")
