### Title:    Multiple Imputation with Bayesian Regularized Regression
### Author:   Kyle M. Lang
### Created:  2014-DEC-12
### Modified: 2016-NOV-05
### Purpose:  The following functions implement MIBEN or MIBL to create multiple
###           imputations within a MICE framework that uses the Bayesian
###           Elastic Net (BEN) or the Bayesian LASSO (BL), respectively, as its
###           elementary imputation method.
### Notes:    The mibrr function implements the imputation models. The R portions
###           of mibrr take care of data pre- and post-processing, while the
###           gibbs sampling and the MCEM optimization of the BEN and BL penalty
###           parameters are done in C++ (see source in runGibbs.cpp). The
###           miben and mibl functions are wrappers that parameterization mibrr
###           as needed to run MIBEN or MIBL, respectively. The ben and bl
###           functions simply fit the Bayesian elastic net and Bayesian LASSO
###           models without any missing data imputation.

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
mibrr <- function(doBl,
                  doImp,
                  data,
                  nImps,
                  targetVars,
                  ignoreVars,
                  iterations,
                  sampleSizes,
                  missCode,
                  returnConvInfo,
                  returnParams,
                  verbose,
                  seed,
                  control)
{
    if(!is.null(seed)) set.seed(seed)
    
    ## Check the user inputs and resolve a set of target variables:
    checkInputs()
    
    ## Define some useful constants and counters:
    nTargets <- length(targetVars)
    nPreds   <- ncol(data) - length(ignoreVars) - 1
    nObs     <- nrow(data)
    
    ## Ensure a correct list of control parameters:
    control <- padControlList()
    
    ## Preprocess the raw data:

    ## Save the original variable names:
    rawNames <- colnames(data)

    ## Set aside the ignored variables:
    ignoredColumns <- data[ , ignoreVars]
    
    ## Move target variables to the leftmost columns
    data <- data.frame(data[ , targetVars],
                       data[ , !colnames(data) %in% c(targetVars, ignoreVars)]
                       )

    ## Hack to deal with 1D matrix conversion to vector:
    if(length(targetVars) == 1) colnames(data)[1] <- targetVars
    
    ## Replace missCode entries with NAs
    if(!is.null(missCode)) {
        userMissCode <- TRUE
        data[data == missCode] <- NA
    } else {
        userMissCode <- FALSE
    }
    
    ## Create a list of missing elements in each target variable
    ## NOTE: Subtract 1 from each index vector to base indices at 0 for C++
    missList <- lapply(data, FUN = function(x) which(is.na(x)) - 1)

    ## Create a vector of response counts:
    if(length(targetVars) > 1)
        respCounts <- colSums(!is.na(data[ , targetVars]))
    else
        respCounts <- sum(!is.na(data[ , targetVars]))
    
#################################################################################
                                        #print("missList[[1]]")
                                        #print(missList[[1]])
    
                                        #print("missList[[2]]")
                                        #print(missList[[2]])
    
                                        #print("missList[[3]]")
                                        #print(missList[[3]])
    
                                        #print("missList[[4]]")
                                        #print(missList[[4]])
    
                                        #rMat <- is.na(data)
#################################################################################
    
    if(control$fimlStarts) {
        ## Singly impute any missing data on auxiliaries:
        imputeCovs()

        ## Mean-center the data using FIML means as centers:
        scaleDataWithFiml()
    } else {
        ## Start the missing values with (temporary) single imputations:
        simpleImpute()
        
        ## Mean-center the data        
        scaleData()
    }
    
    ## Initialize starting values for the Gibbs sampled parameters.
    ## Important to call this before the NAs are replaced with missCode.
    paramStarts <- initializeParams(data     = data, 
                                    nTargets = nTargets,
                                    doBl     = doBl,
                                    control  = control)
    lambdaMat   <- paramStarts$lambda
    betaStarts  <- paramStarts$beta
    tauStarts   <- paramStarts$tau
    sigmaStarts <- paramStarts$sigma

    ## Fill remaining missing data with an integer code:
    if(control$fimlStarts) applyMissCode()
    
#################################################################################
                                        #data[rMat] <- -1e5
#################################################################################
    
    ## Estimate the MIBEN/MIBL model:
    gibbsOut <-
        runGibbs(data            = as.matrix(data),
                 dataScales      = dataScales,
                 nTargets        = nTargets,
                 missList        = missList[c(1 : nTargets)],
                 respCounts      = respCounts,
                 lambda1Starts   = lambdaMat[ , 1], 
                 lambda2Starts   = lambdaMat[ , 2],     # Ignored for BL
                 sigmaStarts     = sigmaStarts,
                 tauStarts       = tauStarts,
                 betaStarts      = betaStarts,
                 nApproxIters    = iterations[1],
                 nTuneIters      = iterations[2],
                 nApproxBurn     = control$approxBurn,
                 nApproxGibbs    = sampleSizes[1],
                 nTuneBurn       = control$tuneBurn,
                 nTuneGibbs      = sampleSizes[2],
                 nPostBurn       = control$postBurn,
                 nPostGibbs      = sampleSizes[3],
                 emConvTol       = control$mcemEpsilon, # Ignored for BL
                 lambdaWindow    = control$smoothingWindow,
                 verbose         = verbose,
                 doBl            = doBl,
                 doImputation    = doImp,
                 adaptScales     = control$adaptScales,
                 simpleIntercept = control$simpleIntercept,
                 twoPhaseOpt     = control$twoPhaseOpt) # Ignored for BL
   
    names(gibbsOut) <- targetVars

    ## Uncenter the data:
    scaleData(revert = TRUE)

    ## Replace missing values:
    if(doImp) {
        missFill <- ifelse(userMissCode, userMissCode, NA)
        for(v in colnames(data)) {
            ## Rebase indices as per R's preference :
            missList[[v]] <- missList[[v]] + 1
            data[missList[[v]], v] <- missFill
        }
    }
    
    ## Compute the potential scale reduction factors (R-Hats) for the posterior
    ## imputation model parameters:
    rHatList <- lapply(X           = c(1 : nTargets),
                       FUN         = checkGibbsConv,
                       gibbsStates = gibbsOut,
                       returnRHats = TRUE,
                       targetNames = targetVars,
                       critVal     = control$convThresh)
    
    ## Draw imputations from the convergent posterior predictive distribution:
    if(doImp) {
        outImps <- getImputedData(gibbsState  = gibbsOut,
                                  nImps       = nImps,
                                  data        = data,
                                  targetMeans = dataMeans[targetVars],
                                  targetVars  = targetVars,
                                  missList    = missList)
    }
    
    ## Aggregate and return the requested output:
    outList <- list()
    if(doImp) outList$imps <- outImps
    
    if(returnConvInfo) {
        outList$rHats <- rHatList
        lamHistList   <- list()
        for(j in 1 : length(gibbsOut)) {
            lamHistList[[j]] <- gibbsOut[[j]]$lambdaHistory
        }
        outList$lambdaHistory <- lamHistList
    }
    
    totalIters <- sum(iterations)
    if(returnParams) {
        betaList <- tauList <- sigmaList <- lamList <- list()
        for(j in 1 : length(gibbsOut)) {
            betaList[[j]]  <- gibbsOut[[j]]$beta
            tauList[[j]]   <- gibbsOut[[j]]$tau
            sigmaList[[j]] <- gibbsOut[[j]]$sigma
            lamList[[j]]   <- gibbsOut[[j]]$lambdaHistory[totalIters, ]
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
miben <- function(data,
                  nImps          = 100,
                  targetVars     = NULL,
                  ignoreVars     = NULL,
                  iterations     = c(100, 10),
                  sampleSizes    = c(25, 250, 500),
                  missCode       = NULL,
                  returnConvInfo = TRUE,
                  returnParams   = FALSE,
                  verbose        = TRUE,
                  seed           = NULL,
                  control        = list()
                  )
{
    mibrr(doBl           = FALSE,
          doImp          = TRUE,
          data           = data,
          nImps          = nImps,
          targetVars     = targetVars,
          ignoreVars     = ignoreVars,
          iterations     = iterations,
          sampleSizes    = sampleSizes,
          missCode       = missCode,
          returnConvInfo = returnConvInfo,
          returnParams   = returnParams,
          verbose        = verbose,
          seed           = seed,
          control        = control)
}# END miben()


### Specify a wrapper function to implement Multiple Imputation with the
### Bayesian Lasso (MIBL):
mibl <- function(data,
                 nImps          = 100,
                 targetVars     = NULL,
                 ignoreVars     = NULL,
                 iterations     = c(100, 10),
                 sampleSizes    = c(25, 250, 500),
                 missCode       = NULL,
                 returnConvInfo = TRUE,
                 returnParams   = FALSE,
                 verbose        = TRUE,
                 seed           = NULL,
                 control        = list()
                 )
{
    mibrr(doBl           = TRUE,
          doImp          = TRUE,
          data           = data,
          nImps          = nImps,
          targetVars     = targetVars,
          ignoreVars     = ignoreVars,
          iterations     = iterations,
          sampleSizes    = sampleSizes,
          missCode       = missCode,
          returnConvInfo = returnConvInfo,
          returnParams   = returnParams,
          verbose        = verbose,
          seed           = seed,
          control        = control)
}# END mibl()


### Specify a wrapper function to fit the Bayesian Elastic Net (BEN):
ben <- function(data,
                y              = NULL,
                X              = NULL,
                iterations     = c(100, 10),
                sampleSizes    = c(25, 250, 500),
                returnConvInfo = TRUE,
                verbose        = TRUE,
                seed           = NULL,
                control        = list()
                )
{
    mibrr(doBl           = FALSE,
          doImp          = FALSE,
          data           = data,
          targetVars     = y,
          ignoreVars     = setdiff(colnames(data), c(y, X)),
          iterations     = iterations,
          sampleSizes    = sampleSizes,
          missCode       = NULL,
          returnConvInfo = returnConvInfo,
          returnParams   = TRUE,
          verbose        = verbose,
          seed           = seed,
          control        = control)
}# END ben()


### Specify a wrapper function to fit the Bayesian Elastic Net (BEN):
bl <- function(data,
               y              = NULL,
               X              = NULL,
               iterations     = c(100, 10),
               sampleSizes    = c(25, 250, 500),
               returnConvInfo = TRUE,
               verbose        = TRUE,
               seed           = NULL,
               control        = list()
               )
{
    mibrr(doBl           = TRUE,
          doImp          = FALSE,
          data           = data,
          targetVars     = y,
          ignoreVars     = setdiff(colnames(data), c(y, X)),
          iterations     = iterations,
          sampleSizes    = sampleSizes,
          missCode       = NULL,
          returnConvInfo = returnConvInfo,
          returnParams   = TRUE,
          verbose        = verbose,
          seed           = seed,
          control        = control)
}# END bl()
