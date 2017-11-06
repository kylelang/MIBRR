### Title:    Test MIBRR Package
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2017-OCT-27
### Purpose:  Script to help test the MIBRR package

rm(list = ls(all = TRUE))

library(mitools)
library(psych)
library(mibrr)

###---------------------------------------------------------------------------###

### Prepare Data for Testing ###

data(bfi)
tmp <- na.omit(bfi)

ed.d           <- model.matrix(~factor(tmp$education))[ , -1]
colnames(ed.d) <-
    c("finish_hs", "some_college", "college_grad", "graduate_degree")

male            <- tmp$gender
male[male == 2] <- 0

cn   <- setdiff(colnames(bfi), c("gender", "education"))
bfi2 <- data.frame(tmp[ , cn], male, ed.d)

rownames(bfi2) <- NULL

targets  <- list(mar = paste0(c("A", "E"), rep(c(1 : 5), each = 2)),
                 mcar = NA,
                 mnar = NA)
pm       <- list(mar = 0.3)
snr      <- list(mar = 5)
marPreds <- c("age",
              "male",
              "finish_hs",
              "some_college",
              "college_grad",
              "graduate_degree")

cn <- c(targets$mar, marPreds)

dat1 <- bfi2[sample(c(1 : nrow(bfi2)), 500), cn]

dat2 <- mibrr::imposeMissing(data    = dat1,
                             targets = targets,
                             preds   = marPreds,
                             pm      = pm,
                             snr     = snr)$data 

###---------------------------------------------------------------------------###

### Check Missing Data Indexing ###

missIndices <- lapply(dat2, function(x) which(is.na(x)) - 1)
scales      <- unlist(lapply(dat2, sd, na.rm = TRUE))
respCounts  <- colSums(!is.na(dat2))

dat3              <- dat2
dat3[is.na(dat2)] <- -9999

for(v in 1 : ncol(dat3)) {
    poiOut <- printObsIndices(data        = as.matrix(dat3),
                              scales      = scales,
                              missIndices = missIndices,
                              respCounts  = respCounts,
                              noMiss      = FALSE,
                              targetIndex = v - 1)

    pmiOut <- printMissIndices(data        = as.matrix(dat3),
                               scales      = scales,
                               missIndices = missIndices,
                               respCounts  = respCounts,
                               noMiss      = FALSE,
                               targetIndex = v - 1)

    cat(paste0("Checking column number ", v, "\n"))
    
    ## Any overlap?
    t1 <- intersect(pmiOut, poiOut)
    if(length(t1) > 0) stop("Miss and obs indices are not disjoint")
    
    ## All rows maintained?
    outInds <- c(pmiOut, poiOut)
    t2 <- setdiff(c(1 : nrow(dat3)), outInds + 1)
    if(length(t2) > 0) stop("Some rows dropped")
    
    ## In == Out?
    t3 <- setdiff(pmiOut, missIndices[[v]])
    if(length(t3) > 0) stop("Miss indices broken")
}# END for(v in 1 : ncol(dat3))

###---------------------------------------------------------------------------###

### Check MIBEN and MIBL ###

nImps <- 25
nReps <- 40

keys <- list(agree = c("-A1", "A2", "A3", "A4", "A5"),
             extra = c("-E1", "-E2", "E3", "E4", "E5")
             )

mibenRes <- miblRes <- miceRes <- list() 
for(rp in 1 : nReps) {
    cat(paste0("Doing rep ", rp, "\n\n"))
    
    dat1 <- bfi2[sample(c(1 : nrow(bfi2)), 500), cn]
    dat2 <- mibrr::imposeMissing(data    = dat1,
                                 targets = targets,
                                 preds   = marPreds,
                                 pm      = pm,
                                 snr     = snr)$data

    cat("Doing MIBEN...")
    mibenOut <- miben(data           = dat2,
                      nImps          = nImps,
                      targetVars     = targets$mar,
                      ignoreVars     = NULL,
                      iterations     = c(50, 10),
                      returnConvInfo = FALSE,
                      returnParams   = FALSE,
                      verbose        = FALSE)
    cat("Done\n")

    cat("Doing MIBL...")
    miblOut <- mibl(data           = dat2,
                    nImps          = nImps,
                    targetVars     = targets$mar,
                    ignoreVars     = NULL,
                    iterations     = c(50, 10),
                    returnConvInfo = FALSE,
                    returnParams   = FALSE,
                    verbose        = FALSE)
    cat("Done\n")

    cat("Doing MICE...")
    miceOut <-
        mice(dat2, m = nImps, maxit = 10, method = "norm", printFlag = FALSE)
    cat("Done\n")

    cat("Fitting models...")
    mibenList <- miblList <- miceList <- list()
    for(m in 1 : nImps) {
        ## MIBEN estimates:
        scores         <- scoreItems(keys  = keys,
                                     items = mibenOut$imps[[m]])$scores
        mibenList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                            mA = mean(scores[ , "agree"]),
                            mE = mean(scores[ , "extra"])
                            )
        ## MIBL estimates:
        scores        <- scoreItems(keys  = keys,
                                    items = miblOut$imps[[m]])$scores
        miblList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                           mA = mean(scores[ , "agree"]),
                           mE = mean(scores[ , "extra"])
                           )
        ## MICE estimates:
        scores        <- scoreItems(keys  = keys,
                                    items = complete(miceOut, m))$scores
        miceList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                           mA = mean(scores[ , "agree"]),
                           mE = mean(scores[ , "extra"])
                           )
    }

    mibenRes[[rp]] <- colMeans(do.call(rbind, mibenList))
    miblRes[[rp]]  <- colMeans(do.call(rbind, miblList))
    miceRes[[rp]]  <- colMeans(do.call(rbind, miceList))
    cat("Done\n\n")
}

mibenFrame <- do.call(rbind, mibenRes)
miblFrame  <- do.call(rbind, miblRes)
miceFrame  <- do.call(rbind, miceRes)

mibenFrame <- rbind(firstRun$miben, mibenFrame)
miblFrame  <- rbind(firstRun$mibl, miblFrame)
miceFrame  <- rbind(firstRun$mice, miceFrame)

simOut <- list(miben = mibenFrame,
               mibl  = miblFrame,
               mice  = miceFrame)

saveRDS(simOut, "miniSimOut.rds")

## Complete data result:
scores  <- scoreItems(keys = keys, items = bfi2)$scores
compRes <- c(r  = cor(scores[ , 1], scores[ , 2]),
             mA = mean(scores[ , "agree"]),
             mE = mean(scores[ , "extra"])
             )

100 * (colMeans(mibenFrame) - compRes) / compRes
100 * (colMeans(miblFrame) - compRes) / compRes
100 * (colMeans(miceFrame) - compRes) / compRes

apply(mibenFrame, 2, sd)
apply(miblFrame, 2, sd)
apply(miceFrame, 2, sd)


###---------------------------------------------------------------------------###

### Test Random Variate Samplers ###

## MVN sampler:
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

## Inverse gamma sampler:
gamShape <- 10
gamScale <- 10

out2.1 <- mibrr::drawInvGamma(nObs, gamShape, gamScale)
out2.2 <- rinvgamma(nObs, gamShape, gamScale)

plot(density(out2.1), col = "red")
lines(density(out2.2), col = "blue")

## Inverse gaussian sampler:
igMu = 1
igLam = 2

out3.1 <- mibrr::drawInvGauss(nObs, igMu, igLam)
out3.2 <- rinvgauss(nObs, igMu, igLam)

plot(density(out3.1), col = "red")
lines(density(out3.2), col = "blue")

## Incomplte gamma calculation:
incGamShape <- 10
incGamCut <- 5

out4.1 <- mibrr::calcIncGamma(incGamShape, incGamCut, FALSE)
out4.2 <- pgamma(q = incGamCut,
                 shape = incGamShape,
                 lower = FALSE) * gamma(incGamShape)

out4.1 - out4.2

###---------------------------------------------------------------------------###

### Test BEN and BL ###

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

###---------------------------------------------------------------------------###

### Check Documentation Examples ###

## Datasets:
data(mibrrExampleData)
data(predictData)

## MIBEN:
mibenOut <- miben(data       = mibrrExampleData,
                  nImps      = 100,
                  iterations = c(30, 10),
                  targetVars = c("y", paste0("x", c(1 : 3))),
                  ignoreVars = "idNum")

## MIBL:
miblOut <- mibl(data       = mibrrExampleData,
                nImps      = 100,
                iterations = c(50, 10),
                targetVars = c("y", paste0("x", c(1 : 3))),
                ignoreVars = "idNum")

## BEN:
benOut <- ben(data       = mibrrExampleData,
              y          = "y",
              X          = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
              iterations = c(30, 10)
              )

## BL:
blOut <- bl(data       = mibrrExampleData,
            y          = "y",
            X          = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
            iterations = c(50, 10)
            )

## predictMibrr:
benOut <- ben(data       = predictData$train,
              y          = "agree",
              X          = setdiff(colnames(predictData$train), "agree"),
              iterations = c(30, 10)
              )
benPred <- predictMibrr(object = benOut, newData = predictData$test)

mibenOut <- miben(data         = predictData$incomplete,
                  nImps        = 100,
                  iterations   = c(30, 10),
                  returnParams = TRUE)
mibenPred <- predictMibrr(object = mibenOut, newData = predictData$test)

## Info Functions:
mibrrL()
mibrrW()
