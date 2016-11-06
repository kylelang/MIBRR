### Title:    MIBEN Simulation Functions
### Author:   Kyle M. Lang
### Created:  2016-MAY-05
### Modified: 2016-MAY-06

                                        #rm(list = ls(all = TRUE))

prepData <- function() {
    data("MASchools")
    
    ## Clean the raw MA Schools data:
    dat1              <- MASchools[ , -c(1, 2)]
    expNames          <- colnames(dat1)[grep("exp", colnames(dat1))]
    dat1[ , expNames] <- dat1[ , expNames] / 100
    
    targets  <- c("score8", "exptot", "salary", "stratio")
    nTargets <- length(targets)
    auxVars  <- setdiff(colnames(dat1),
                        c(targets, "score4", "expbil", "expocc", "english"))
    
    dat2 <- na.omit(dat1[ , c(targets, auxVars)])
    
    rawCor        <- cor(dat2)
    rawCov        <- cov(dat2)
    rawSds        <- sqrt(diag(rawCov))
    auxCor        <- rawCor[auxVars, auxVars]
    aux2TargetCor <- rawCor[auxVars, targets]

    targetMu <- colMeans(dat2[ , targets])
    auxMeans <-colMeans(dat2[ , auxVars])
       
    ## Get Sigma for "small data" version:
    smallSigma <-
        cov(dat2[ , c("score8", "exptot", "stratio", "salary", "income")],
            use = "pairwise")
    smallMu <-
        colMeans(dat2[ , c("score8", "exptot", "stratio", "salary", "income")],
                 na.rm = TRUE)
    
    ## Return the necessary data moments:
    list(
        smallSigma        = smallSigma,
        smallMu           = smallMu,
        targetSigma       = rawCov[targets, targets],
        targetMu          = targetMu,
        auxScaleMean      = mean(rawSds[auxVars]),
        auxScaleSd        = sd(rawSds[auxVars]),
        auxCorMean        = mean(auxCor),
        auxCorSd          = sd(auxCor),
        aux2TargetCorMean = mean(aux2TargetCor),
        aux2TargetCorSd   = sd(aux2TargetCor),
        auxMeanMean       = mean(auxMeans),
        auxMeanSd         = sd(auxMeans)
    )
}# END prepData()


startRng <- function(rp, parms) {
    if(!rp %in% .lec.GetStreams()) {
        .lec.CreateStream(c(1 : parms$maxStreams))
        .lec.CurrentStream(rp)
    } else {
        .lec.ResetStartStream(rp)
        .lec.CurrentStream(rp)
    }
}


getSigma <- function(nAux, moments) {
    ## Unpack parameters:
    auxCorMean        <- moments$auxCorMean
    auxCorSd          <- moments$auxCorSd
    aux2TargetCorMean <- moments$aux2TargetCorMean
    aux2TargetCorSd   <- moments$aux2TargetCorSd
    auxScaleMean      <- moments$auxScaleMean
    auxScaleSd        <- moments$auxScaleSd
    targetSigma       <- moments$targetSigma
    targetSds         <- sqrt(diag(targetSigma))
    nTargets          <- length(moments$targetMu)
    nVars             <- nTargets + nAux
    
    badSigma <- TRUE
    while(badSigma) {
        ## Simulate a matrix of auxiliary variable correlations:
        tmpMat <- matrix(0, nAux, nAux)
        for(i in 1 : nAux)
            for(j in i : nAux)
                tmpMat[i, j] <- rnorm(1, auxCorMean, auxCorSd)
        
        auxCor <- tmpMat + t(tmpMat)
        diag(auxCor) <- 1.0
        
        ## Simulate a matrix of auxiliary to target variable correlations:
        aux2TargetCor <-
            matrix(rnorm(nAux * nTargets, aux2TargetCorMean, aux2TargetCorSd),
                   nAux,
                   nTargets)
        
        ## Simulate a vector of auxiliary variable scales:
        auxSdShape <- auxScaleMean^2 / auxScaleSd^2
        auxSdScale <- auxScaleSd^2 / auxScaleMean
        auxSds <- rgamma(nAux, shape = auxSdShape, scale = auxSdScale)

        ## Construct weight matrices to convert the correlation matrix to a
        ## covariance matrix:
        m1 <- matrix(c(targetSds, auxSds), nVars, nVars)
        m2 <- matrix(c(targetSds, auxSds), nVars, nVars, byrow = TRUE)

        ## Populate the full correlation matrix:
        sigma <- matrix(NA, nVars, nVars)
        sigma[nTargets + 1 : nAux, nTargets + 1 : nAux] <- auxCor
        sigma[1 : nTargets, nTargets + 1 : nAux] <- t(aux2TargetCor)
        sigma[nTargets + 1 : nAux, 1 : nTargets] <- aux2TargetCor
        sigma <- sigma * m1 * m2 # Convert correlations to variances/covariances
        sigma[1 : nTargets, 1 : nTargets] <- targetSigma
        
        ## Ensure that the full covariance matrix is positive definite:
        eVals <- eigen(sigma)$values
        if(any(eVals < 0)) {
            correctSigma <- nearPD(sigma)
            sigma <- as.matrix(correctSigma$mat)
            badSigma <- !correctSigma$converged
        } else {
            badSigma <- FALSE
        }
    }# CLOSE while(badSigma)          
    rownames(sigma) <- colnames(sigma) <-
        c(names(targetMu), paste0("z", c(1 : nAux)))
    sigma
}# END getSigma()


getMu <- function(nAux, moments) {
    auxMean <- moments$auxMeanMean
    auxSd <- moments$auxMeanSd
    targetMu <- moments$targetMu
    
    tmp <- rnorm(nAux, auxMean, auxSd)
    names(tmp) <- paste0("z", c(1 : nAux))
    c(targetMu, tmp)
}# END getMu()


simulateData <- function(parms) {
    if(parms$bigData) {
        sigma   <- getSigma(parms$nAux, parms$moments)
        mu      <- getMu(parms$nAux, parms$moments)
    } else {
        sigma <- parms$moments$smallSigma
        mu    <- parms$moments$smallMu
    }
    simData <- as.data.frame(rmvnorm(parms$nObs, mu, sigma))
    simData
}# END simData()


imposeMissing <- function(rawData, parms)
{
    incompVars <- parms$incompVars
    auxVars    <- parms$auxVars
    pm         <- parms$pm
    marType    <- parms$marType
    
    if(length(auxVars) > 1) auxVar <- rowMeans(rawData[ , auxVars])
    else                    auxVar <- rawData[ , auxVars]
    
    pVec <- pnorm(auxVar, mean(auxVar), sd(auxVar))
    
    for(v in 1 : length(incompVars)) {
        if(marType[v] == "tails") {
            rVec <- pVec < (pm / 2) | pVec > (1 - (pm/2))
        } else if(marType[v] == "center") {
            rVec <- pVec > (0.5 - (pm / 2)) & pVec < (0.5 + (pm / 2))
        } else if(marType[v] == "lower") {
            rVec <- pVec < pm
        } else if(marType[v] == "upper") {
            rVec <- pVec > (1 - pm)
        } else {
            stop("Please provide a valid 'marType'")
        }
        rawData[rVec, incompVars[v]] <- NA
    }
    rawData
}# END imposeMissing()


getIters <- function(parms)
{
    switch(parms$pm * 10,
           list(
               approxIters = 100,
               tuneIters   = 10,
               approxN     = 25,
               smooth      = 10
           ),
           list(
               approxIters = 150,
               tuneIters   = 20,
               approxN     = 50,
               smooth      = 15
           ),
           list(
               approxIters = 200,
               tuneIters   = 30,
               approxN     = 100,
               smooth      = 20
           ),
           list(
               approxIters = 300,
               tuneIters   = 40,
               approxN     = 150,
               smooth      = 30
           )
           )
}# END getIters()


doRep <- function(rp, parms) {
    cat(paste0("Doing replication ", rp, "\n"))
    
    ## Initialize the RNG:
    startRng(rp = rp, parms = parms)
    
    testForm <- parms$testForm
    resDir   <- parms$resDir
    
    for(pm in parms$pmVec) {# LOOP over proportions missing
        parms$pm <- pm
        resList <- rHatList <- lambdaList <- list()
        
### Data Generation ###
        simData  <- simulateData(parms = parms)
        
        ## Select auxiliary variables:
        if(expNum == 1) parms$auxVars <- "income"
        else parms$auxVars <- paste0("z", sample(c(1 : parms$nAux), 10))
        
        missData <- imposeMissing(rawData = simData, parms = parms)
        
### Missing Data Imputation ###
        iterList <- getIters(parms = parms)
        mibenOut <- try(
            miben(rawData         = missData,
                  targetVars      = parms$incompVars,
                  nImps           = parms$nImps,
                  verboseIters    = FALSE,
                  verboseErrors   = FALSE,
                  mcemApproxIters = iterList$approxIters,
                  mcemTuneIters   = iterList$tuneIters,
                  mcemApproxN     = iterList$approxN,
                  control         = list(smoothingWindow = iterList$smooth)
                  )
        )
        
        miblOut <- try(
            mibl(rawData         = missData,
                 targetVars      = parms$incompVars,
                 nImps           = parms$nImps,
                 verboseIters    = FALSE,
                 verboseErrors   = FALSE,
                 mcemApproxIters = iterList$approxIters,
                 mcemTuneIters   = iterList$tuneIters,
                 mcemApproxN     = iterList$approxN,
                 control         = list(smoothingWindow = iterList$smooth)
                 )
        )
        
        miceOut <- try(
            mice(data      = missData,
                 m         = parms$nImps,
                 maxit     = 5,
                 printFlag = FALSE)
        )
        
### Model Fitting ###
        ## Fit complete data models:
        resList$comp <- lm(testForm, data = simData)
        
        ## Fit MIBEN models:
        if(class(mibenOut) != "try-error") {
            fitList <-
                lapply(mibenOut$imps, FUN = function(x) lm(testForm, data = x))
            resList$miben <- MIcombine(fitList)
            rHatList$miben <- mibenOut$rHats
            lambdaList$miben <- mibenOut$lambdaHistory
        } else {
            resList$miben <- mibenOut
        }
        
        ## Fit MIBL models:
        if(class(miblOut) != "try-error") {
            fitList <-
                lapply(miblOut$imps, FUN = function(x) lm(testForm, data = x))
            resList$mibl <- MIcombine(fitList)
            rHatList$mibl <- miblOut$rHats
            lambdaList$mibl <- miblOut$lambdaHistory
        } else {
            resList$mibl <- miblOut
        }
        
        ## Fit MICE models:
        if(class(miceOut) != "try-error") {
            impList <- list()
            for(m in 1 : parms$nImps) impList[[m]] <- complete(miceOut, m)
            
            fitList <- lapply(impList, FUN = function(x) lm(testForm, data = x))
            resList$mice <- MIcombine(fitList)
        } else {
            resList$mice <- miceOut
        }

        saveRDS(iterList, paste0(resDir,
                                 "mibenSimIters",
                                 "_exp", expNum,
                                 "_n",   parms$nObs,
                                 "_v",   parms$nAux,
                                 "_pm",  100*parms$pm,
                                 "_rep", rp,
                                 ".rds")
                )
        
        saveRDS(resList, paste0(resDir,
                                "mibenSimOut",
                                "_exp", expNum,
                                "_n",   parms$nObs,
                                "_v",   parms$nAux,
                                "_pm",  100*parms$pm,
                                "_rep", rp,
                                ".rds")
                )
        
        saveRDS(rHatList, paste0(resDir,
                                 "mibenSimRHats",
                                 "_exp", expNum,
                                 "_n",   parms$nObs,
                                 "_v",   parms$nAux,
                                 "_pm",  100*parms$pm,
                                 "_rep", rp,
                                 ".rds")
                )
        
        saveRDS(lambdaList, paste0(resDir,
                                   "mibenSimLambda",
                                   "_exp", expNum,
                                   "_n",   parms$nObs,
                                   "_v",   parms$nAux,
                                   "_pm",  100*parms$pm,
                                   "_rep", rp,
                                   ".rds")
                )
    }# END for(pm in parms$pmVec)
    
    list()# Return nothing
}# END doRep()


### Summary Functions for the Simulation Results ###

getResults <- function(nReps, parms)
{
    nTargets <- length(parms$incompVars)
    nPreds   <- parms$nAux + nTargets - 1
    
    failCount <- 0
    for(rp in 1 : nReps) {
        tmp <- readRDS(paste0(parms$resDir,
                              "mibenSimOut",
                              "_exp", parms$expNum,
                              "_n",   parms$nObs,
                              "_v",   parms$nAux,
                              "_pm",  100*parms$pm,
                              "_rep", rp,
                              ".rds")
                       )
        
        tmp2 <- readRDS(paste0(parms$resDir,
                               "mibenSimRHats",
                               "_exp", parms$expNum,
                               "_n",   parms$nObs,
                               "_v",   parms$nAux,
                               "_pm",  100*parms$pm,
                               "_rep", rp,
                               ".rds")
                        )
        
        if(rp == 1) {
            tmpList <- list()
            for(i in names(tmp))
                tmpList[[i]] <-
                    matrix(NA, nReps, nTargets,
                           dimnames = list(NULL, names(coef(tmp$comp)))
                           )
            coefList <- seList <- tmpList
            
            tmpList <- rHatList <- list()
            for(i in names(tmp2)) {
                rHatNames <-
                    unlist(lapply(parms$incompVars,
                                  FUN = function(x, y) paste(x, y, sep = ":"),
                                  y = names(unlist(tmp2$miben[[1]]))
                                  )
                           )
                rHatList[[i]] <- matrix(NA, nReps, 2 * (nPreds + 1) * nTargets,
                                        dimnames = list(NULL, rHatNames)
                                        )
                tmpList[[i]] <- matrix(NA, nReps, nTargets)
            }
            fmiList <- dfList <- tmpList
        }
        
        if(!is.null(tmp) & class(tmp) != "character") {
            for(i in names(tmp)) {
                coefList[[i]][rp, ]  <- coef(tmp[[i]])            
                seList[[i]][rp, ]  <- sqrt(diag(vcov(tmp[[i]])))
            }
            for(i in names(tmp2)) {
                fmiList[[i]][rp, ] <- tmp[[i]]$missinfo            
                dfList[[i]][rp, ] <- tmp[[i]]$df           
                rHatList[[i]][rp, ] <- unlist(tmp2[[i]])
            }
        } else {
            failCount <- failCount + 1
        }
    }
    list(nFails = failCount,
         rHats  = rHatList,
         coef   = coefList,
         se     = seList,
         fmi    = fmiList,
         df     = dfList)
}# END getResults()


calcBias <- function(coefList, compName = "comp")
{
    prbList <- sbList <- list()
    for(i in setdiff(names(coefList), compName)) {
        coefs <- coefList[[i]]
        means <- colMeans(coefs, na.rm = TRUE)
        compMeans <- colMeans(coefList[[compName]], na.rm = TRUE)
        
        prbList[[i]] <- 100 * (means - compMeans) / abs(compMeans)
        sbList[[i]]  <- (means - compMeans) / apply(coefs, 2, sd)
    }
    list(prb = prbList,
         sb  = sbList)
}


tCrit <- function(x, alpha) qt(p = 1 - (alpha / 2), df = x)


calcCiStats <- function(coefList,
                        seList,
                        dfList,
                        alpha = 0.05,
                        compName = "comp")
{
    cicList <- ciwList <- list()
    for(i in names(dfList)) {
        coefs  <- coefList[[i]]
        ses    <- seList[[i]]
        dfs    <- dfList[[i]]
        tCrits <- apply(dfs, c(1, 2), tCrit, alpha = alpha)
        
        moes <- tCrits * ses
        lbs  <- coefs - moes
        ubs  <- coefs + moes
        
        tmp <- coefList[[compName]] < ubs & coefList[[compName]] > lbs
        cicList[[i]] <- colMeans(tmp)

        ciwList[[i]] <- colMeans(ubs - lbs)
    }
    list(cic = cicList,
         ciw = ciwList)
}



testFun <- function(rp, parms) {
    print(paste0("Doing replication ", rp))
    
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
    
    testOut <- ben(rawData         = dat1,
                   y               = "y",
                   X               = paste0("x", c(1 : nPreds)),
                   mcemApproxIters = parms$approxIters,
                   mcemApproxN     = parms$approxN,
                   mcemTuneIters   = parms$tuneIters,
                   mcemTuneN       = parms$tuneN,
                   mcemPostN       = parms$postN,
                   verbose         = parms$verbose,
                   control         =
                       list(center          = parms$center,
                            scale           = parms$scale,
                            adaptScales     = FALSE,
                            regIntercept    = FALSE,
                            useClassic      = parms$useClassic,
                            simpleIntercept = parms$simpleInt)
                   )
    
    testOut2 <- bl(rawData         = dat1,
                   y               = "y",
                   X               = paste0("x", c(1 : nPreds)),
                   mcemApproxIters = parms$approxIters,
                   mcemApproxN     = parms$approxN,
                   mcemTuneIters   = parms$tuneIters,
                   mcemTuneN       = parms$tuneN,
                   mcemPostN       = parms$postN,
                   verbose         = parms$verbose,
                   control         =
                       list(center          = parms$center,
                            scale           = parms$scale,
                            adaptScales     = FALSE,
                            regIntercept    = FALSE,
                            useClassic      = parms$useClassic,
                            simpleIntercept = parms$simpleInt)
                   )

    testForm <-
        as.formula(paste0("y ~ ", paste0("x", c(1 : nPreds), collapse = " + ")))
    lmOut <- lm(testForm, data = dat1)
    
    testDat <- simulateData(nObs   = nObs,
                            nPreds = nPreds,
                            r2     = r2,
                            collin = collin,
                            beta   = beta,
                            means  = means,
                            scales = scales)

    if(parms$center) {
        testDat2 <- data.frame(scale(testDat, scale = FALSE))
        extra <- mean(dat1$y)
    } else {
        testDat2 <- testDat
        extra <- 0.0
    }
    if(parms$scale) {
        tmpSds    <- apply(testDat2, 2, sd)
        tmpSds[1] <- 1.0
        testDat2  <- data.frame(
            scale(testDat2, center = FALSE, scale = tmpSds)
        )
    }
    
    predOut1 <- predictMibrr(object = testOut,
                             newData = as.matrix(testDat2[ , -1])
                             ) + extra
    predOut2 <- predictMibrr(object = testOut2,
                             newData = as.matrix(testDat2[ , -1])
                             ) + extra
    predOut3 <- predict(lmOut, newdata = testDat)
    
    mseVec <- c(
        mean((predOut1 - dat1$y)^2),
        mean((predOut2 - dat1$y)^2),
        mean((predOut3 - dat1$y)^2)
    )
    mseVec
}# END testFun()
