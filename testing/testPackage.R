### Title:    Test MIBRR Package
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2019-JAN-23

rm(list = ls(all = TRUE))

library(mitools)
library(psych)
library(MIBRR)
library(devtools)
library(parallel)
library(MCMCpack)
library(statmod)
library(HyperbolicDist)
library(LaplacesDemon)

install_github("kylelang/SURF/source/SURF")
library(SURF)

MIBRR:::testMissIndex()
MIBRR:::testSamplers()

###--------------------------------------------------------------------------###

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
cn       <- c(targets$mar, marPreds)

dat1 <- bfi2[sample(c(1 : nrow(bfi2)), 500), cn]
dat2 <- imposeMissData(data    = dat1,
                       targets = targets,
                       preds   = marPreds,
                       pm      = pm,
                       snr     = snr)$data 

###--------------------------------------------------------------------------###

### Check Missing Data Indexing ###

missIndices <- lapply(dat2, function(x) which(is.na(x)) - 1)
scales      <- unlist(lapply(dat2, sd, na.rm = TRUE))
respCounts  <- colSums(!is.na(dat2))

dat3              <- dat2
dat3[is.na(dat2)] <- -9999

for(v in 1 : ncol(dat3)) {
    poiOut <- MIBRR:::printObsIndices(data        = as.matrix(dat3),
                                      scales      = scales,
                                      missIndices = missIndices,
                                      respCounts  = respCounts,
                                      noMiss      = FALSE,
                                      targetIndex = v - 1)

    pmiOut <- MIBRR:::printMissIndices(data        = as.matrix(dat3),
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

###--------------------------------------------------------------------------###

### Test MIBEN and MIBL ###

testFun <- function(rp, data, env) {
    cn       <- env$cn
    targets  <- env$targets
    marPreds <- env$marPreds
    pm       <- env$pm
    snr      <- env$snr
    nImps    <- env$nImps
    keys     <- env$keys
    
    dat1 <- data[sample(c(1 : nrow(data)), 500), cn]
    dat2 <- imposeMissData(data    = dat1,
                           targets = targets,
                           preds   = marPreds,
                           pm      = pm,
                           snr     = snr)$data

    mibenOut <- miben(data       = dat2,
                      targetVars = targets$mar,
                      ignoreVars = NULL,
                      iterations = c(50, 10),
                      verbose    = FALSE)
    mibenImps <- getImpData(mibenOut, nImps)
    
    miblOut <- mibl(data       = dat2,
                    targetVars = targets$mar,
                    ignoreVars = NULL,
                    iterations = c(50, 10),
                    verbose    = FALSE)
    miblImps <- getImpData(miblOut, nImps)
    
    miceOut <-
        mice(dat2, m = nImps, maxit = 10, method = "norm", printFlag = FALSE)

    mibenList <- miblList <- miceList <- list()
    for(m in 1 : nImps) {
        ## MIBEN estimates:
        scores         <- scoreItems(keys  = keys,
                                     items = mibenImps[[m]])$scores
        mibenList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                            mA = mean(scores[ , "agree"]),
                            mE = mean(scores[ , "extra"])
                            )
        ## MIBL estimates:
        scores        <- scoreItems(keys  = keys,
                                    items = miblImps[[m]])$scores
        miblList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                           mA = mean(scores[ , "agree"]),
                           mE = mean(scores[ , "extra"])
                           )
        ## MICE estimates:
        scores        <- scoreItems(keys  = keys,
                                    items = mice::complete(miceOut, m))$scores
        miceList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                           mA = mean(scores[ , "agree"]),
                           mE = mean(scores[ , "extra"])
                           )
    }

    list(miben = colMeans(do.call(rbind, mibenList)),
         mibl  = colMeans(do.call(rbind, miblList)),
         mice  = colMeans(do.call(rbind, miceList))
         )
} # END testFun()

nReps <- 4
nImps <- 10
keys  <- list(agree = c("-A1", "A2", "A3", "A4", "A5"),
              extra = c("-E1", "-E2", "E3", "E4", "E5")
              )

simOut <- mclapply(X        = c(1 : nReps),
                   FUN      = testFun,
                   data     = bfi2,
                   env      = parent.frame(),
                   mc.cores = 4)

tmp <- do.call(rbind, lapply(simOut, unlist))

mibenFrame <- tmp[ , grep("miben", colnames(tmp))]
miblFrame  <- tmp[ , grep("mibl", colnames(tmp))]
miceFrame  <- tmp[ , grep("mice", colnames(tmp))]

## Complete data result: 
scores  <- scoreItems(keys = keys, items = bfi2)$scores
compRes <- c(r  = cor(scores[ , 1], scores[ , 2]),
             mA = mean(scores[ , "agree"]),
             mE = mean(scores[ , "extra"])
             )

## Percent Relative Bias:
100 * (colMeans(mibenFrame) - compRes) / compRes
100 * (colMeans(miblFrame) - compRes) / compRes
100 * (colMeans(miceFrame) - compRes) / compRes

# Monte Carlo SD:
apply(mibenFrame, 2, sd)
apply(miblFrame, 2, sd)
apply(miceFrame, 2, sd)

###--------------------------------------------------------------------------###

### Test Vanilla MI ###

testFun <- function(rp, data, env) {
    cn       <- env$cn
    targets  <- env$targets
    marPreds <- env$marPreds
    pm       <- env$pm
    snr      <- env$snr
    nImps    <- env$nImps
    keys     <- env$keys
    
    dat1 <- data[sample(c(1 : nrow(data)), 500), cn]
    dat2 <- imposeMissData(data    = dat1,
                           targets = targets,
                           preds   = marPreds,
                           pm      = pm,
                           snr     = snr)$data
    
    vanOut <- vanilla(data       = dat2,
                      targetVars = targets$mar,
                      ignoreVars = NULL,
                      verbose    = FALSE)
    
    vanImps <- getImpData(vanOut, nImps)
    
    miceOut <-
        mice(dat2, m = nImps, maxit = 10, method = "norm", printFlag = FALSE)

    vanList <- miceList <- list()
    for(m in 1 : nImps) {
        ## Vanilla estimates:
        scores         <- scoreItems(keys  = keys,
                                     items = vanImps[[m]])$scores
        vanList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                          mA = mean(scores[ , "agree"]),
                          mE = mean(scores[ , "extra"])
                          )

        ## MICE estimates:
        scores        <- scoreItems(keys  = keys,
                                    items = mice::complete(miceOut, m))$scores
        miceList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                           mA = mean(scores[ , "agree"]),
                           mE = mean(scores[ , "extra"])
                           )
    }

    list(van  = colMeans(do.call(rbind, vanList)),
         mice = colMeans(do.call(rbind, miceList))
         )
} # END testFun()

nReps <- 20
nImps <- 10
keys  <- list(agree = c("-A1", "A2", "A3", "A4", "A5"),
              extra = c("-E1", "-E2", "E3", "E4", "E5")
              )

simOut <- mclapply(X        = c(1 : nReps),
                   FUN      = testFun,
                   data     = bfi2,
                   env      = parent.frame(),
                   mc.cores = 4)

tmp <- do.call(rbind, lapply(simOut, unlist))

vanFrame  <- tmp[ , grep("van", colnames(tmp))]
miceFrame <- tmp[ , grep("mice", colnames(tmp))]

## Complete data result: 
scores  <- scoreItems(keys = keys, items = bfi2)$scores
compRes <- c(r  = cor(scores[ , 1], scores[ , 2]),
             mA = mean(scores[ , "agree"]),
             mE = mean(scores[ , "extra"])
             )

## Percent Relative Bias:
100 * (colMeans(vanFrame) - compRes) / compRes
100 * (colMeans(miceFrame) - compRes) / compRes

# Monte Carlo SD:
apply(vanFrame, 2, sd)
apply(miceFrame, 2, sd)

###--------------------------------------------------------------------------###

### Test Random Variate Samplers ###

## MVN sampler:
nObs           <- 500000
mvnMu          <- rep(10, 3)
mvnSigma       <- matrix(5, 3, 3)
diag(mvnSigma) <- 20

out1.1 <- MIBRR:::drawMvn(nObs, mvnMu, mvnSigma, seed = 235711)
out1.2 <- rmvnorm(nObs, mvnMu, mvnSigma)

par(mfrow = c(1, 3))
for(i in 1 : length(mvnMu)) {
    plot(density(out1.1[ , i]), col = "red")
    lines(density(out1.2[ , i]), col = "blue")
}

ks.test(x = out1.1[ , 3], y = out1.2[ , 3])

## Inverse gamma sampler:
gamShape <- 10
gamScale <- 10

out2.1 <- MIBRR:::drawInvGamma(nObs, gamShape, gamScale, 235711)
out2.2 <- rinvgamma(nObs, gamShape, gamScale)

plot(density(out2.1), col = "red")
lines(density(out2.2), col = "blue")

ks.test(x = out2.1, y = out2.2)

## Inverse gaussian sampler:
igMu  <- 1
igLam <- 2

out3.1 <- MIBRR:::drawInvGauss(nObs, igMu, igLam, 235711)
out3.2 <- rinvgauss(nObs, igMu, igLam)

plot(density(out3.1), col = "red")
lines(density(out3.2), col = "blue")

ks.test(x = out3.1, y = out3.2)

## GIG sampler:
gigLam <- 1
gigChi <- 2
gigPsi <- 2

out4.1 <- MIBRR:::drawGig(nObs, gigLam, gigChi, gigPsi, 235711)
out4.2 <- rgig(nObs, c(gigLam, gigChi, gigPsi))

plot(density(out4.1), col = "red")
lines(density(out4.2), col = "blue")

ks.test(x = out4.1, y = out4.2)

## Scaled inverse chi-squared sampler:
df    <- 100
scale <- 10

out5.1 <- MIBRR:::drawScaledInvChiSq(nObs, df, scale, 235711)
out5.2 <- rinvchisq(nObs, df, scale)

plot(density(out5.1), col = "red")
lines(density(out5.2), col = "blue")

ks.test(x = out4.1, y = out4.2)

## Incomplete gamma calculation:
incGamShape <- 10
incGamCut   <- 5

out6.1 <- MIBRR:::calcIncGamma(incGamShape, incGamCut, FALSE)
out6.2 <- pgamma(q     = incGamCut,
                 shape = incGamShape,
                 lower = FALSE) * gamma(incGamShape)

all.equal(out6.1, out6.2)

###--------------------------------------------------------------------------###

### Test BEN and BL ###

testFun <- function(rp, parms) {
    nObs   <- parms$nObs
    nPreds <- parms$nPreds
    r2     <- parms$r2
    collin <- parms$collin
    beta   <- parms$beta
    means  <- parms$means
    scales <- parms$scales
    ipp    <- parms$ipp
    lamSq  <- parms$reliability
    nIvs   <- ifelse(ipp > 1, nPreds * ipp, nPreds)
    
    dat1 <- simRegData(nObs            = nObs,
                       nPreds          = nPreds,
                       r2              = r2,
                       collin          = collin,
                       beta            = beta,
                       means           = means,
                       scales          = scales,
                       itemsPerPred    = ipp,
                       predReliability = lamSq)
    
    enOut <- try(
        ben(data    = dat1,
            y       = "y",
            X       = paste0("x", c(1 : nIvs)),
            verbose = parms$verbose,
            control = list(center          = parms$center,
                           scale           = parms$scale,
                           adaptScales     = parms$adaptScales,
                           simpleIntercept = parms$simpleInt,
                           optCheckKkt     = parms$checkKkt)
            )
    )
    
    blOut <- try(
        bl(data    = dat1,
           y       = "y",
           X       = paste0("x", c(1 : nIvs)),
           verbose = parms$verbose,
           control = list(center          = parms$center,
                          scale           = parms$scale,
                          adaptScales     = parms$adaptScales,
                          simpleIntercept = parms$simpleInt)
           )
    )

    form1 <-
        as.formula(paste0("y ~ ", paste0("x", c(1 : nIvs), collapse = " + ")))
    lmOut <- try(lm(form1, data = dat1))
    
    mseMat <- matrix(NA, parms$nTests, 3)
    for(i in 1 : parms$nTests) {
        testDat <- simRegData(nObs            = nObs,
                              nPreds          = nPreds,
                              r2              = r2,
                              collin          = collin,
                              beta            = beta,
                              means           = means,
                              scales          = scales,
                              itemsPerPred    = ipp,
                              predReliability = lamSq)
        
        if(class(enOut) != "try-error") {
            enPred       <- postPredict(mibrrFit = enOut, newData = testDat)$y
            mseMat[i, 1] <- mean((enPred - testDat$y)^2)
        } else {
            mseMat[i, 1] <- NA
        }
        
        if(class(blOut) != "try-error") {
            blPred       <- postPredict(mibrrFit = blOut, newData = testDat)$y
            mseMat[i, 2] <- mean((blPred - testDat$y)^2)
        } else {
            mseMat[i, 2] <- NA
        }
        
        if(class(lmOut) != "try-error") {
            lmPred       <- predict(lmOut, newdata = testDat)
            mseMat[i, 3] <- mean((lmPred - testDat$y)^2)
        } else {
            mseMat[i, 3] <- NA
        }
    }
       
    outMat           <- rbind(colMeans(mseMat), apply(mseMat, 2, var))
    colnames(outMat) <- c("BEN", "BL", "MLR")
    rownames(outMat) <- c("MSE", "var(MSE)")

    errorList <- list()
    if(class(enOut) == "try-error") errorList$ben <- enOut
    if(class(blOut) == "try-error") errorList$bl  <- blOut
    if(class(lmOut) == "try-error") errorList$lm  <- lmOut
    
    list(results = outMat, errors = errorList)
}# END testFun()

## Parameterize mini-simulation:
alpha  <- 0.5
nPreds <- 8

parms             <- list()
parms$nObs        <- 100
parms$nPreds      <- nPreds
parms$r2          <- 0.5
parms$collin      <- 0.1
parms$beta        <- matrix(c(alpha, runif(nPreds, 0.3, 0.6)))
parms$means       <- runif(nPreds, 0, 1)
parms$scales      <- rep(1, nPreds)
parms$center      <- TRUE
parms$scale       <- TRUE
parms$simpleInt   <- FALSE
parms$verbose     <- FALSE
parms$adaptScales <- TRUE
parms$nTests      <- 100
parms$checkKkt    <- FALSE
parms$ipp         <- 10
parms$reliability <- 0.8

## Run mini-simulation in parallel:
nReps   <- 8
mseList <- mclapply(X        = c(1 : nReps),
                    FUN      = testFun,
                    parms    = parms,
                    mc.cores = 4)

mseList

###--------------------------------------------------------------------------###

### Check Documentation Examples ###

## MIBEN ##
?miben

data(mibrrExampleData)

## MCEM estimation:
mibenOut <- miben(data       = mibrrExampleData,
                  iterations = c(30, 10),
                  targetVars = c("y", paste0("x", c(1 : 3))),
                  ignoreVars = "idNum")

## Fully Bayesian estimation:
mibenOut <- miben(data         = mibrrExampleData,
                  targetVars   = c("y", paste0("x", c(1 : 3))),
                  ignoreVars   = "idNum",
                  sampleSizes  = c(500, 500),
                  doMcem       = FALSE,
                  lam1PriorPar = c(1.0, 0.1),
                  lam2PriorPar = c(1.0, 0.1)
                  )

## MIBL ##
?mibl

data(mibrrExampleData)

## MCEM estimation:
miblOut <- mibl(data       = mibrrExampleData,
                iterations = c(50, 10),
                targetVars = c("y", paste0("x", c(1 : 3))),
                ignoreVars = "idNum")

## Fully Bayesian estimation:
miblOut <- mibl(data         = mibrrExampleData,
                targetVars   = c("y", paste0("x", c(1 : 3))),
                ignoreVars   = "idNum",
                sampleSizes  = c(500, 500),
                doMcem       = FALSE,
                lam1PriorPar = c(1.0, 0.1)
                )

## BEN ##
?ben

data(mibrrExampleData)

## MCEM estimation:
benOut <- ben(data       = mibrrExampleData,
              y          = "y",
              X          = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
              iterations = c(30, 10)
              )

## Fully Bayesian estimation:
benOut <- ben(data         = mibrrExampleData,
              y            = "y",
              X            = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
              doMcem       = FALSE,
              sampleSizes  = c(500, 500),
              lam1PriorPar = c(1.0, 0.1),
              lam2PriorPar = c(1.0, 0.1)
              )
     
## BL ##
?bl

data(mibrrExampleData)

## MCEM estimation:
blOut <- bl(data       = mibrrExampleData,
            y          = "y",
            X          = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
            iterations = c(50, 10)
            )

## Fully Bayesian estimation:
blOut <- bl(data         = mibrrExampleData,
            y            = "y",
            X            = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
            doMcem       = FALSE,
            sampleSizes  = c(500, 500),
            lam1PriorPar = c(1.0, 0.1)
            )
     
## postPredict ##
?postPredict

data(predictData)

## Fit a Bayesian elastic net model:
benOut <- ben(data       = predictData$train,
              y          = "agree",
              X          = setdiff(colnames(predictData$train), "agree"),
              iterations = c(30, 10)
              )

## Generate 10 posterior predictions for 'y':
benPred <- postPredict(mibrrFit = benOut,
                       newData  = predictData$test,
                       nDraws   = 10)

## Generate posterior predictions for 'y' using MAP scores (the default):
benPred <- postPredict(mibrrFit = benOut,
                       newData  = predictData$test,
                       nDraws   = 0)

## Generate posterior predictions for 'y' using EAP scores:
benPred <- postPredict(mibrrFit = benOut,
                       newData  = predictData$test,
                       nDraws   = -1)

## Fit chained Bayesian elastic net models to support MI:
mibenOut <- miben(data = predictData$incomplete, iterations = c(30, 10))

## Generate posterior predictions from elementary imputation models:
mibenPred <- postPredict(mibrrFit = mibenOut, newData = predictData$test)

## getImpData ##
?getImpData

data(mibrrExampleData)

## Fit a Bayesian elastic net model:
mibenOut <- miben(data         = mibrrExampleData,
                  targetVars   = c("y", paste0("x", c(1 : 3))),
                  ignoreVars   = "idNum",
                  sampleSizes  = c(500, 500),
                  doMcem       = FALSE,
                  lam1PriorPar = c(1.0, 0.1),
                  lam2PriorPar = c(1.0, 0.1)
                  )

## Generate 25 imputed datasets:
mibenImps <- getImpData(mibrrFit = mibenOut, nImps = 25)

## getParams ##
?getParams

data(mibrrExampleData)

## Fit a Bayesian elastic net model:
benOut <- ben(data         = mibrrExampleData,
              y            = "y",
              X            = paste0("x", c(1 : 3)),
              sampleSizes  = c(500, 500),
              doMcem       = FALSE,
              lam1PriorPar = c(1.0, 0.1),
              lam2PriorPar = c(1.0, 0.1)
              )

## Extract posterior parameter samples:
benPars <- getParams(mibrrFit = benOut, target = "y")

## getField ##
?getField

data(mibrrExampleData)

## Fit a Bayesian elastic net model:
benOut <- ben(data         = mibrrExampleData,
              y            = "y",
              X            = paste0("x", c(1 : 3)),
              sampleSizes  = c(500, 500),
              doMcem       = FALSE,
              lam1PriorPar = c(1.0, 0.1),
              lam2PriorPar = c(1.0, 0.1)
              )

## Extract the potential scale reduction factors:
benRHats <- getField(mibrrFit = benOut, what = "rHats")

## Info Functions ##
?mibrrL
mibrrL()

?mibrrW
mibrrW()
