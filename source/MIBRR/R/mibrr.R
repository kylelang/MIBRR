### Title:    Multiple Imputation with Bayesian Regularized Regression
### Author:   Kyle M. Lang
### Created:  2014-DEC-12
### Modified: 2017-NOV-17
### Purpose:  The following functions implement MIBEN or MIBL to create multiple
###           imputations within a MICE framework that uses the Bayesian
###           Elastic Net (BEN) or the Bayesian LASSO (BL), respectively, as its
###           elementary imputation method.
### Notes:    - The mibrr function implements the imputation models.
###           - The R layer of mibrr take care of data pre- and post-processing
###             as well as the MCEM optimization of the BEN and BL penalty.
###           - The gibbs sampling is done in C++ (see source in runGibbs.cpp).
###           - The miben and mibl functions are wrappers that parameterization
###             mibrr as needed to run MIBEN or MIBL, respectively.
###           - The ben and bl functions simply fit the Bayesian elastic net and
###             Bayesian LASSO models to (possibly) incomplete data without
###             returning any missing data imputations.

##--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------##
##  Copyright (C) 2017 Kyle M. Lang <k.m.lang@uvt.nl>                          ##  
##                                                                             ##
##  This file is part of MIBRR.                                                ##
##                                                                             ##
##  This program is free software: you can redistribute it and/or modify it    ##
##  under the terms of the GNU General Public License as published by the      ##
##  Free Software Foundation, either version 3 of the License, or (at you      ##
##  option) any later version.                                                 ##
##                                                                             ##
##  This program is distributed in the hope that it will be useful, but        ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of                 ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General   ##
##  Public License for more details.                                           ##
##                                                                             ##
##  You should have received a copy of the GNU General Public License along    ##
##  with this program. If not, see <http://www.gnu.org/licenses/>.             ##
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
    
    ## Save the original variable names:
    rawNames <- colnames(data)

    ## Set aside the ignored variables:
    ignoredColumns <- as.data.frame(data[ , ignoreVars])
    if(length(ignoreVars) == 1) colnames(ignoredColumns) <- ignoreVars
    
    ## Move target variables to the leftmost columns
    data <- data.frame(data[ , targetVars],
                       data[ , !colnames(data) %in% c(targetVars, ignoreVars)]
                       )

    ## Hack to deal with 1D matrix conversion to vector:
    if(length(targetVars) == 1) colnames(data)[1] <- targetVars
    
    ## Replace missCode entries with NAs
    if(!is.null(missCode)) {
        userMissCode           <- TRUE
        data[data == missCode] <- NA
    } else {
        userMissCode <- FALSE
    }
    
    ## Create a list of missing elements in each target variable
    ## NOTE: Subtract 1 from each index vector to base indices at 0 for C++
    missList <- lapply(data, FUN = function(x) which(is.na(x)) - 1)
      
    ## Create a vector of response counts:
    respCounts <- colSums(!is.na(data))
    noMiss     <- all(respCounts == nObs)

    dataMeans <- dataScales <- rep(NA, ncol(data))
    if(noMiss) {
        ## Mean-center the data (and compute initial data scales):
        scaleData()
    } else if(control$fimlStarts) {
        ## Singly impute any missing data on auxiliaries:
        imputeCovs()
        ## Mean-center the data using FIML means as centers:
        scaleDataWithFiml()
    } else {
        ## Start the missing values with (temporary) single imputations:
        simpleImpute()
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
    if(control$fimlStarts & !noMiss) applyMissCode()

    ## Calculate the total number of MCEM iterations:
    totalIters <- sum(iterations)
    
    ## Create a container for Lambda's iteration history:
    lambdaHistory <- lapply(targetVars,
                            function(x) {
                                tmp <- matrix(NA, totalIters - 1, 2)
                                colnames(tmp) <- c("lambda1", "lambda2")
                                tmp
                            }
                            )
    names(lambdaHistory) <- targetVars
    
    for(i in 1 : totalIters) {
        ## Print status update:
        if(verbose) {
            if(i == 1) {
                cat("\nBeginning MCEM 'Approximation' phase\n")
                mcemStage <- "'Warm-Up'"
            }
            if(i == (iterations[1] + 1)){
                cat("\nBeginning MCEM 'Tuning' phase\n")
                mcemStage <- "'Tuning'"
            }
            if(i == totalIters) cat("\nSampling from the stationary posterior\n")
        }
        
        ## What Gibbs sample sizes should we use?:
        if     (i <= iterations[1]) sams <- sampleSizes[[1]]
        else if(i < totalIters    ) sams <- sampleSizes[[2]]
        else                        sams <- sampleSizes[[3]]
        
        if(verbose) cat("\n") # Beautify output
        
        ## Estimate the MIBEN/MIBL model:
        gibbsOut <-
            runGibbs(data            = as.matrix(data),
                     dataScales      = dataScales,
                     nTargets        = nTargets,
                     missList        = missList[c(1 : nTargets)],
                     respCounts      = respCounts[c(1 : nTargets)],
                     lambda1         = lambdaMat[ , 1], 
                     lambda2         = lambdaMat[ , 2], # Ignored for BL
                     sigmaStarts     = sigmaStarts,
                     tauStarts       = tauStarts,
                     betaStarts      = betaStarts,
                     burnSams        = sams[1],
                     totalSams       = sum(sams),
                     verbose         = verbose,
                     doBl            = doBl,
                     adaptScales     = control$adaptScales,
                     simpleIntercept = control$simpleIntercept,
                     noMiss          = noMiss)

        if(i < totalIters) {
            if(verbose) {
                check <- i > iterations[1]
                cat(paste0("Doing MCEM ",
                           ifelse(check, "'Tuning'", "'Approximation'"),
                           " iteration ",
                           ifelse(check, i - iterations[1], i),
                           " of ",
                           iterations[as.numeric(check) + 1],
                           "\n")
                    )
            }
            
            ## Conduct the MCEM update of the lambdas:
            optOut <- optimizeLambda(lambdaMat    = lambdaMat,
                                     gibbsState   = gibbsOut,
                                     doBl         = doBl,
                                     returnACov   = FALSE, # Do we ever care?
                                     controlParms = list(
                                         method      = control$optMethod,
                                         boundLambda = control$optBoundLambda,
                                         showWarns   = verbose,
                                         traceLevel  = control$optTraceLevel,
                                         checkKkt    = control$optCheckKkt
                                     )
                                     )
            
            if(doBl) {
                lambdaMat           <- cbind(unlist(optOut), NA)
                colnames(lambdaMat) <- c("lambda", "dummy")
            } else {
                if(control$optCheckKkt) {
                    optCols <- 4
                    optLabs <- c("kkt1", "kkt2", "lambda1", "lambda2")
                } else {
                    optCols <- 2
                    optLabs <- c("lambda1", "lambda2")
                }
                
                ## Organize the optimization output:
                optMat <- matrix(unlist(optOut),
                                 ncol     = optCols,
                                 byrow    = TRUE,
                                 dimnames = list(NULL, optLabs)
                                 )
                lambdaMat <- optMat[ , grep("lambda", colnames(optMat))]

                ## Hack for models with single DV:
                if(length(targetVars) == 1) lambdaMat <- matrix(lambdaMat, 1, 2)
                
                ## Check optimality conditions:
                if(control$optCheckKkt) {
                    kkt1Flag <- optMat[ , "kkt1"] == 0
                    kkt2Flag <- optMat[ , "kkt2"] == 0
                    
                    if(any(kkt1Flag))
                        stop("First KKT optimality condition not satisfied when optimizing Lambda")
                    
                    if(any(kkt2Flag))
                        stop("Second KKT optimality condition not satisfied when optimizing Lambda")
                }
            }

            for(j in 1 : nTargets) {
                betaStarts[ , j] <- colMeans(gibbsOut[[j]]$beta[ , -1])
                sigmaStarts[j]   <- mean(gibbsOut[[j]]$sigma)
                tauStarts[ , j]  <- colMeans(gibbsOut[[j]]$tau)
                
                lambdaHistory[[j]][i, ] <- lambdaMat[j, ]

                if(i == iterations[1]) {
                    smoothRange    <- (i - control$smoothingWindow + 1) : i
                    lambdaMat[j, ] <- colMeans(lambdaHistory[[j]][smoothRange, ])
                }
            }
        }
    }# END for(i in 1 : totalIters)
    
    ## Give some nicer names:
    names(gibbsOut) <- rownames(lambdaMat) <- targetVars

    ## Uncenter the data:
    scaleData(revert = TRUE)

    ## Replace missing values:
    if(doImp) {
        missFill <- ifelse(userMissCode, userMissCode, NA)
        for(v in colnames(data)) {
            ## Rebase indices as per R's preference :
            missList[[v]]          <- missList[[v]] + 1
            data[missList[[v]], v] <- missFill
        }
    }
    
    ## Compute the potential scale reduction factors (R-Hats) for the posterior
    ## imputation model parameters:
    rHatList <- lapply(X           = targetVars,
                       FUN         = checkGibbsConv,
                       gibbsStates = gibbsOut,
                       returnRHats = TRUE,
                       critVal     = control$convThresh)
    
    ## Draw imputations from the convergent posterior predictive distribution:
    if(doImp) outImps <- getImputedData(gibbsState  = gibbsOut,
                                        nImps       = nImps,
                                        data        = data,
                                        targetMeans = dataMeans[targetVars],
                                        targetVars  = targetVars,
                                        missList    = missList)
    
    ## Provide some pretty names for the output objects:
    nameOutput()
    
    ## Aggregate and return the requested output:
    outList <- list()
    if(doImp) outList$imps <- outImps
    
    if(returnConvInfo) {
        outList$rHats         <- rHatList
        outList$lambdaHistory <- lambdaHistory
    }
    
    if(returnParams) {
        paramList <- list()
        for(j in targetVars) {
            paramList[[j]]$beta   <- gibbsOut[[j]]$beta
            paramList[[j]]$tau    <- gibbsOut[[j]]$tau
            paramList[[j]]$sigma  <- gibbsOut[[j]]$sigma
            paramList[[j]]$lambda <- lambdaMat[j, ]
        }
        outList$params <- paramList
    }
    outList
}# END mibrr()


### Specify a wrapper function to implement Multiple Imputation with the
### Bayesian Elastic Net (MIBEN):
miben <- function(data,
                  nImps,
                  targetVars     = NULL,
                  ignoreVars     = NULL,
                  iterations     = c(100, 10),
                  sampleSizes    = list(rep(25, 2), rep(250, 2), rep(500, 2)),
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
                 nImps,
                 targetVars     = NULL,
                 ignoreVars     = NULL,
                 iterations     = c(100, 10),
                 sampleSizes    = list(rep(25, 2), rep(250, 2), rep(500, 2)),
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
                y,
                X              = NULL,
                iterations     = c(100, 10),
                sampleSizes    = list(rep(25, 2), rep(250, 2), rep(500, 2)),
                missCode       = NULL,
                returnConvInfo = TRUE,
                verbose        = TRUE,
                seed           = NULL,
                control        = list()
                )
{
    if(length(y) > 1) stop("Only one outcome variable is allowed.")
    
    mibrr(doBl           = FALSE,
          doImp          = FALSE,
          data           = data,
          targetVars     = y,
          ignoreVars     = setdiff(colnames(data), c(y, X)),
          iterations     = iterations,
          sampleSizes    = sampleSizes,
          missCode       = missCode,
          returnConvInfo = returnConvInfo,
          returnParams   = TRUE,
          verbose        = verbose,
          seed           = seed,
          control        = control)
}# END ben()


### Specify a wrapper function to fit the Bayesian LASSO (BL):
bl <- function(data,
               y,
               X              = NULL,
               iterations     = c(100, 10),
               sampleSizes    = list(rep(25, 2), rep(250, 2), rep(500, 2)),
               missCode       = NULL,
               returnConvInfo = TRUE,
               verbose        = TRUE,
               seed           = NULL,
               control        = list()
               )
{
    if(length(y) > 1) stop("Only one outcome variable is allowed.")
    
    mibrr(doBl           = TRUE,
          doImp          = FALSE,
          data           = data,
          targetVars     = y,
          ignoreVars     = setdiff(colnames(data), c(y, X)),
          iterations     = iterations,
          sampleSizes    = sampleSizes,
          missCode       = missCode,
          returnConvInfo = returnConvInfo,
          returnParams   = TRUE,
          verbose        = verbose,
          seed           = seed,
          control        = control)
}# END bl()
