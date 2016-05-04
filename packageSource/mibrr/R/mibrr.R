### Title:    Multiple Imputation with Bayesian Regularized Regression
### Author:   Kyle M. Lang
### Created:  2014-DEC-12
### Modified: 2016-MAY-04
### Purpose:  The following functions implements MIBEN or MIBL to create multiple
###           imputations within a MICE framework that uses the Bayesian
###           Elastic Net (BEN) or the Bayesian LASSO (BL), respectively, as its
###           elementary imputation method.
### Notes:    The mibrr function implements the imputation models. The R portions
###           of mibrr take care of data pre- and post-processing, while the
###           gibbs sampling and the MCEM optimization of the BEN and BL penalty
###           parameters are done in C++ (see source in runGibbs.cpp). The
###           miben and mibl functions are wrappers that parameterization mibrr
###           as needed to run MIBEN or MIBL, respectively.

##--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------##
##  Copyright (C) 2016 Kyle M. Lang <kyle.lang@ttu.edu>                        ##  
##                                                                             ##
##  This file is part of mibrr.                                                ##
##                                                                             ##
##  This program is free software: you can redistribute it and/or modify it    ##
##  under the terms of the GNU Lesser General Public License as published by   ##
##  the Free Software Foundation, either version 3 of the License, or          ##
##  (at you option) any later version.                                         ##
##                                                                             ##
##  This program is distributed in the hope that it will be useful, but        ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY ##
##  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public    ##
##  License for more details.                                                  ##
##                                                                             ##
##  You should have received a copy of the GNU Lesser General Public License   ##
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.      ##
##-----------------------------------------------------------------------------##

## Specify the main computational function of the mibrr package:
mibrr <- function(doMiben,
                  rawData,
                  nImps             = 100,
                  targetVariables   = NULL,
                  ignoreVariables   = NULL,
                  nEmApproxIters    = 100,
                  nEmTuneIters      = 10,
                  nEmApproxBurn     = 25,
                  nEmApproxGibbs    = 25,
                  nEmTuneBurn       = 100,
                  nEmTuneGibbs      = 200,
                  nPosteriorBurn    = NULL,
                  nPosteriorThin    = 5,
                  missCode          = NULL,
                  returnConvInfo    = TRUE,
                  returnModelParams = FALSE,
                  verboseIters      = TRUE,
                  verboseErrors     = TRUE,
                  gibbsControl      = list(),
                  optControl        = list()
                  )
{
    ## Populate the control parameters:
    gibbsControlDefaults =
        list(createRngStream = TRUE,
             streamName = NULL,
             convThresh = 1.1,
             rngSeed = 235711)
    
    if(doMiben) {
        optControlDefaults =
            list(lambda1Starts = NULL,
                 lambda2Starts = NULL,
                 emEpsilon = 1.0e-5,
                 smoothingWindow = 1)
    } else {
        optControlDefaults =
            list(lambdaStarts = NULL,
                 usePCStarts = FALSE,
                 smoothingWindow = 1)
    }
    
    gibbsControl <- lapply(names(gibbsControlDefaults),
                           FUN = padControlList,
                           inList = gibbsControl,
                           defaultList = gibbsControlDefaults)
    names(gibbsControl) <- names(gibbsControlDefaults)
    
    optControl <- lapply(names(optControlDefaults),
                         FUN = padControlList,
                         inList = optControl,
                         defaultList = optControlDefaults)
    names(optControl) <- names(optControlDefaults)
    
    ## Check for target variables. When no targets are given, all incomplete
    ## variables not listed in 'ignoreVariables' are imputed.
    if(is.null(targetVariables)) {
        targetCandidates <-
            colnames(rawData)[!colnames(rawData) %in% ignoreVariables]
        warning("You did not specify any target variables, so I will impute \
the missing data on\nevery variable in 'rawData' that is not listed in \
'ignoreVariables'.\n")        
    } else {
        targetCandidates <- targetVariables
    }
    
    ## Make sure 'rawData' contains missing data that we can find:
    if(is.null(missCode)){
        completeTargets <- colMeans(is.na(rawData[ , targetCandidates])) == 0
        if(all(completeTargets)) {
            stop("Your target variables appear to be fully observed. Did you \
forget to provide a\nvalue for 'missCode'?\n")
        }
    } else {
        rMat <- rawData == missCode
        if(!any(rMat, na.rm = TRUE)) {
            stop(paste0("The value you provided for 'missCode' (i.e., ",
                        missCode,
                        ") does not appear anywhere in 'rawData'.\n",
                        "Are you sure that ",
                        missCode,
                        " encodes your missing data?\n")
                 )
        } else {
            rawData[rMat] <- NA
        }
    }

    ## Select the final set of target variables:
    targetVars <- targetCandidates[!completeTargets]
    if(any(completeTargets)) {
        warning(paste0("The potential target variables {",
                       paste(targetCandidates[completeTargets],
                             collapse = ", "),
                       "} are fully observed.\n",
                       "These items will not be imputed.\n")
                )
    }
    
    ## Define some useful constants and counters:
    nTargets <- length(targetVars)
    nPreds <- ncol(rawData) - length(ignoreVariables) - 1
    nObs <- nrow(rawData)
    
    ## Populate the starting values for Lambda:
    if(doMiben) {
        if(is.null(optControl$lambda1Starts)) {
            lambda1Starts <- rep(0.5, nTargets)
        } else {
            options(warn = -1)# Suppress warnings about recycling elements
            lambda1Starts <- matrix(optControl$lambda1Starts, nTargets, 1)
            options(warn = 0)
        }
        if(is.null(optControl$lambda2Starts)) {
            lambda2Starts <- rep(nPreds / 10, nTargets)
        } else {
            options(warn = -1)# Suppress warnings about recycling elements
            lambda2Starts <- matrix(optControl$lambda2Starts, nTargets, 1)
            options(warn = 0)
        }
        lambdaMat <- cbind(lambda1Starts, lambda2Starts)
    } else {
        if(is.null(optControl$lambdaStarts)) {
            if(optControl$usePCStarts) {
                ## Must call this before the NA's are replaced with missCode:
                lambdaVec <- getLambdaStarts(inData = centeredData,
                                             nTargets = nTargets,
                                             nPreds = nPreds)
            } else {
                lambdaVec <- rep((nPreds / 10), nTargets)
                lambdaMat <- cbind(lambdaVec, 0)
            }
        } else {
            options(warn = -1)
            lambdaVec <- as.vector(matrix(optControl$lambdaStarts, nTargets, 1))
            options(warn = 0)
            lambdaMat <- cbind(lambdaVec, 0)
        }
    }# END if(doMiben)
    
    ## Preprocess the raw data:
    ## 1) Move target variables to the leftmost columns
    ## 2) Temporarily replace missCode entries with NAs
    ## 3) Mean-center the data
    tmpData <- data.frame(rawData[ , targetVars],
                          rawData[ , !colnames(rawData) %in%
                                  c(targetVars, ignoreVariables)]
                          )
    
    if(!is.null(missCode)) {
        tmpData[tmpData == missCode] <- NA
    }
    
    centeredData <- scale(tmpData, scale = FALSE)
    targetMeans <- attr(centeredData, "scaled:center")[1 : nTargets]
    dataScales <- apply(centeredData, 2, FUN = sd, na.rm = TRUE)
     
    ## Initialize starting values for the Gibbs sampled parameters.
    ## Important to call this before the NAs are replaced with missCode.
    paramStarts <- initializeParams(rawData = centeredData,
                                    inLambda = lambdaMat,
                                    nTargets = nTargets,
                                    doMiben = doMiben)
    betaStarts <- paramStarts$beta
    tauStarts <- paramStarts$tau
    sigmaStarts <- paramStarts$sigma
    
    ## Construct an integer-valued missing data code that
    ## does not take legal data values and use it to flag NAs.
    if(is.null(missCode)) {
        if(max(abs(centeredData), na.rm = TRUE) < 1.0) {
            missCode <- -9
        } else {
            codeMag <- floor(log10(max(abs(centeredData), na.rm = TRUE))) + 2
            missCode <- -(10^codeMag - 1)
        }
    }
    centeredData[is.na(centeredData)] <- missCode
    
    gibbsOut <-
        runGibbs(inData                         = centeredData,
                 dataMeans                      = targetMeans,
                 dataScales                     = dataScales,
                 nTargets                       = nTargets,
                 lambda1Starts                  = lambdaMat[ , 1],
                 lambda2Starts                  = lambdaMat[ , 2],
                 sigmaStarts                    = sigmaStarts,
                 tauStarts                      = tauStarts,
                 betaStarts                     = betaStarts,
                 missCode                       = missCode,
                 nEmApproxIters                 = nEmApproxIters,
                 nEmTuneIters                   = nEmTuneIters,
                 emApprox_nBurnIns              = nEmApproxBurn,
                 emApprox_gibbsSampleSize       = nEmApproxGibbs,
                 emTune_nBurnIns                = nEmTuneBurn,
                 emTune_gibbsSampleSize         = nEmTuneGibbs,
                 posteriorGibbs_nBurnIns        =
                     ifelse(is.null(nPosteriorBurn), -1, nPosteriorBurn),
                 posteriorGibbs_gibbsSampleSize = (nPosteriorThin * nImps),
                 emConvTol                      =
                     ifelse(doMiben, optControl$emEpsilon, -1),
                 lambdaWindow                   = optControl$smoothingWindow,
                 verboseIters                   = verboseIters,
                 verboseErrors                  = verboseErrors,
                 doMibl                         = !doMiben)

    names(gibbsOut) <- targetVars
    
    ## Compute the potential scale reduction factors (R-Hats) for the posterior
    ## imputation model parameters:
    rHatList <- lapply(c(1 : nTargets),
                       FUN         = checkGibbsConv,
                       gibbsStates = gibbsOut,
                       returnRHats = TRUE,
                       targetNames = targetVars,
                       critVal     = gibbsControl$convThresh)
    
    ## Draw imputations from the convergent posterior predictive distribution:
    outImps <- getImputedData(gibbsState  = gibbsOut,
                              nImps       = nImps,
                              rawData     = rawData,
                              targetVars  = targetVars,
                              targetMeans = targetMeans)
    
    ## Aggregate and return the requested output:
    outList <- list()
    outList$imps <- outImps
    
    if(returnConvInfo) {
        outList$rHats <- rHatList
        lamHistList <- list()
        for(j in 1 : length(gibbsOut)) {
            lamHistList[[j]] <- gibbsOut[[j]]$lambdaHistory
        }
        outList$lambdaHistory <- lamHistList
    }
    
    totalEmIters <- nEmApproxIters + nEmTuneIters
    if(returnModelParams) {
        betaList <- tauList <- sigmaList <- lamList <- list()
        for(j in 1 : length(gibbsOut)) {
            betaList[[j]]  <- gibbsOut[[j]]$beta
            tauList[[j]]   <- gibbsOut[[j]]$tau
            sigmaList[[j]] <- gibbsOut[[j]]$sigma
            lamList[[j]]   <- gibbsOut[[j]]$lambdaHistory[totalEmIters, ]
        }
        outList$params <- list(
            beta   = betaList,
            tau    = tauList,
            sigma  = sigmaList,
            lambda = lamList)
    }
    outList
}# END mibrr()



### Specify a wrapper function to implement Multiple Imputation with the
### Bayesian Elastic Net (MIBEN):
miben <- function(rawData,
                  nImps             = 100,
                  targetVariables   = NULL,
                  ignoreVariables   = NULL,
                  nEmApproxIters    = 100,
                  nEmTuneIters      = 10,
                  nEmApproxBurn     = 25,
                  nEmApproxGibbs    = 25,
                  nEmTuneBurn       = 100,
                  nEmTuneGibbs      = 200,
                  nPosteriorBurn    = NULL,
                  nPosteriorThin    = 5,
                  missCode          = NULL,
                  returnConvInfo    = TRUE,
                  returnModelParams = FALSE,
                  verboseIters      = TRUE,
                  verboseErrors     = TRUE,
                  gibbsControl      = list(
                      createRngStream = TRUE,
                      streamName      = NULL,
                      convThresh      = 1.1,
                      rngSeed         = 235711
                  ),
                  optControl        = list(
                      lambda1Starts   = NULL,
                      lambda2Starts   = NULL,
                      emEpsilon       = 1.0e-5,
                      smoothingWindow = 1)
                  )
{
    mibrr(doMiben           = TRUE,
          rawData           = rawData,
          nImps             = nImps,
          targetVariables   = targetVariables,
          ignoreVariables   = ignoreVariables,
          nEmApproxIters    = nEmApproxIters,
          nEmTuneIters      = nEmTuneIters,
          nEmApproxBurn     = nEmApproxBurn,
          nEmApproxGibbs    = nEmApproxGibbs,
          nEmTuneBurn       = nEmTuneBurn,
          nEmTuneGibbs      = nEmTuneGibbs,
          nPosteriorBurn    = nPosteriorBurn,
          nPosteriorThin    = nPosteriorThin,
          missCode          = missCode,
          returnConvInfo    = returnConvInfo,
          returnModelParams = returnModelParams,
          verboseIters      = verboseIters,
          verboseErrors     = verboseErrors,
          gibbsControl      = gibbsControl,
          optControl        = optControl)
}# END miben()



### Specify a wrapper function to implement Multiple Imputation with the
### Bayesian Lasso (MIBL):
mibl <- function(rawData,
                 nImps             = 100,
                 targetVariables   = NULL,
                 ignoreVariables   = NULL,
                 nEmApproxIters    = 100,
                 nEmTuneIters      = 10,
                 nEmApproxBurn     = 25,
                 nEmApproxGibbs    = 25,
                 nEmTuneBurn       = 100,
                 nEmTuneGibbs      = 200,
                 nPosteriorBurn    = NULL,
                 nPosteriorThin    = 5,
                 missCode          = NULL,
                 returnConvInfo    = TRUE,
                 returnModelParams = FALSE,
                 verboseIters      = TRUE,
                 verboseErrors     = TRUE,
                 gibbsControl      = list(
                     createRngStream = TRUE,
                     streamName      = NULL,
                     convThresh      = 1.1,
                     rngSeed         = 235711
                 ),
                 optControl        = list(
                     lambdaStarts    = NULL,
                     usePCStarts     = FALSE,
                     smoothingWindow = 1)
                 )
{
    mibrr(doMiben           = FALSE,
          rawData           = rawData,
          nImps             = nImps,
          targetVariables   = targetVariables,
          ignoreVariables   = ignoreVariables,
          nEmApproxIters    = nEmApproxIters,
          nEmTuneIters      = nEmTuneIters,
          nEmApproxBurn     = nEmApproxBurn,
          nEmApproxGibbs    = nEmApproxGibbs,
          nEmTuneBurn       = nEmTuneBurn,
          nEmTuneGibbs      = nEmTuneGibbs,
          nPosteriorBurn    = nPosteriorBurn,
          nPosteriorThin    = nPosteriorThin,
          missCode          = missCode,
          returnConvInfo    = returnConvInfo,
          returnModelParams = returnModelParams,
          verboseIters      = verboseIters,
          verboseErrors     = verboseErrors,
          gibbsControl      = gibbsControl,
          optControl        = optControl)
}# END mibl()
