### Title:    Multiple Imputation with Bayesian Regularized Regression
### Author:   Kyle M. Lang
### Created:  2014-DEC-12
### Modified: 2016-MAY-05
### Purpose:  The following functions implement MIBEN or MIBL to create multiple
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
                  nImps,
                  targetVars,
                  ignoreVars,
                  mcemApproxIters,
                  mcemTuneIters,
                  mcemApproxN,
                  mcemTuneN,
                  mcemPostN,
                  missCode,
                  returnConvInfo,
                  returnParams,
                  verboseIters,
                  verboseErrors,
                  seed,
                  control)
{
    if(!is.null(seed)) set.seed(seed)
    
    ## Check the user inputs and resolve a set of target variables:
    checkInputs()
    
    ## Define some useful constants and counters:
    nTargets <- length(targetVars)
    nPreds   <- ncol(rawData) - length(ignoreVars) - 1
    nObs     <- nrow(rawData)
    
    ## Ensure a correct list of control parameters:
    control <- padControlList()
    
    ## Preprocess the raw data:
    ## 1) Move target variables to the leftmost columns
    ## 2) Temporarily replace missCode entries with NAs
    ## 3) Mean-center the data
    tmpData <- data.frame(rawData[ , targetVars],
                          rawData[ , !colnames(rawData) %in%
                                  c(targetVars, ignoreVars)]
                          )
    
    if(!is.null(missCode)) tmpData[tmpData == missCode] <- NA
    
    centeredData <- scale(tmpData, scale = FALSE)
    targetMeans  <- attr(centeredData, "scaled:center")[1 : nTargets]
    dataScales   <- apply(centeredData, 2, FUN = sd, na.rm = TRUE)
     
    ## Initialize starting values for the Gibbs sampled parameters.
    ## Important to call this before the NAs are replaced with missCode.
    paramStarts <- initializeParams(rawData  = centeredData,
                                    nTargets = nTargets,
                                    doMiben  = doMiben,
                                    control  = control)
    lambdaMat   <- paramStarts$lambda
    betaStarts  <- paramStarts$beta
    tauStarts   <- paramStarts$tau
    sigmaStarts <- paramStarts$sigma

    ## Fill the missing data with an integer code:
    applyMissCode()
    
    ## Estimate the MIBEN/MIBL model:
    gibbsOut <-
        runGibbs(inData           = centeredData,
                 dataMeans        = targetMeans,
                 dataScales       = dataScales,
                 nTargets         = nTargets,
                 lambda1Starts    = lambdaMat[ , 1],
                 lambda2Starts    = lambdaMat[ , 2],
                 sigmaStarts      = sigmaStarts,
                 tauStarts        = tauStarts,
                 betaStarts       = betaStarts,
                 missCode         = missCode,
                 nMcemApproxIters = mcemApproxIters,
                 nMcemTuneIters   = mcemTuneIters,
                 nMcemApproxBurn  = control$mcemApproxBurn,
                 nMcemApproxGibbs = mcemApproxN,
                 nMcemTuneBurn    = control$mcemTuneBurn,
                 nMcemTuneGibbs   = mcemTuneN,
                 nMcemPostBurn    = ifelse(is.null(control$mcemPostBurn), -1,
                     control$mcemPostBurn),
                 nMcemPostGibbs   = mcemPostN,
                 emConvTol        = ifelse(doMiben, control$mcemEpsilon, -1),
                 lambdaWindow     = control$smoothingWindow,
                 verboseIters     = verboseIters,
                 verboseErrors    = verboseErrors,
                 doMibl           = !doMiben)

    names(gibbsOut) <- targetVars
    
    ## Compute the potential scale reduction factors (R-Hats) for the posterior
    ## imputation model parameters:
    rHatList <- lapply(c(1 : nTargets),
                       FUN         = checkGibbsConv,
                       gibbsStates = gibbsOut,
                       returnRHats = TRUE,
                       targetNames = targetVars,
                       critVal     = control$convThresh)
    
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
    
    totalMcemIters <- mcemApproxIters + mcemTuneIters
    if(returnParams) {
        betaList <- tauList <- sigmaList <- lamList <- list()
        for(j in 1 : length(gibbsOut)) {
            betaList[[j]]  <- gibbsOut[[j]]$beta
            tauList[[j]]   <- gibbsOut[[j]]$tau
            sigmaList[[j]] <- gibbsOut[[j]]$sigma
            lamList[[j]]   <- gibbsOut[[j]]$lambdaHistory[totalMcemIters, ]
        }
        outList$params <- list(
            beta   = betaList,
            tau    = tauList,
            sigma  = sigmaList,
            lambda = lamList
        )
    }
    outList
}# END mibrr()


### Specify a wrapper function to implement Multiple Imputation with the
### Bayesian Elastic Net (MIBEN):
miben <- function(rawData,
                  nImps           = 100,
                  targetVars      = NULL,
                  ignoreVars      = NULL,
                  mcemApproxIters = 100,
                  mcemTuneIters   = 10,
                  mcemApproxN     = 25,
                  mcemTuneN       = 250,
                  mcemPostN       = 500,
                  missCode        = NULL,
                  returnConvInfo  = TRUE,
                  returnParams    = FALSE,
                  verboseIters    = TRUE,
                  verboseErrors   = TRUE,
                  seed            = NULL,
                  control         = list()
                  )
{
    mibrr(doMiben         = TRUE,
          rawData         = rawData,
          nImps           = nImps,
          targetVars      = targetVars,
          ignoreVars      = ignoreVars,
          mcemApproxIters = mcemApproxIters,
          mcemTuneIters   = mcemTuneIters,
          mcemApproxN     = mcemApproxN,
          mcemTuneN       = mcemTuneN,
          mcemPostN       = mcemPostN,
          missCode        = missCode,
          returnConvInfo  = returnConvInfo,
          returnParams    = returnParams,
          verboseIters    = verboseIters,
          verboseErrors   = verboseErrors,
          seed            = seed,
          control         = control)
}# END miben()


### Specify a wrapper function to implement Multiple Imputation with the
### Bayesian Lasso (MIBL):
mibl <- function(rawData,
                 nImps           = 100,
                 targetVars      = NULL,
                 ignoreVars      = NULL,
                 mcemApproxIters = 100,
                 mcemTuneIters   = 10,
                 mcemApproxN     = 25,
                 mcemTuneN       = 250,
                 mcemPostN       = 500,
                 missCode        = NULL,
                 returnConvInfo  = TRUE,
                 returnParams    = FALSE,
                 verboseIters    = TRUE,
                 verboseErrors   = TRUE,
                 seed            = NULL,
                 control         = list()
                 )
{
    mibrr(doMiben         = FALSE,
          rawData         = rawData,
          nImps           = nImps,
          targetVars      = targetVars,
          ignoreVars      = ignoreVars,
          mcemApproxIters = mcemApproxIters,
          mcemTuneIters   = mcemTuneIters,
          mcemApproxN     = mcemApproxN,
          mcemTuneN       = mcemTuneN,
          mcemPostN       = mcemPostN,
          missCode        = missCode,
          returnConvInfo  = returnConvInfo,
          returnParams    = returnParams,
          verboseIters    = verboseIters,
          verboseErrors   = verboseErrors,
          seed            = seed,
          control         = control)
}# END mibl()
