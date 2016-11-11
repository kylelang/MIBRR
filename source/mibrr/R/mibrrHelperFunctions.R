### Title:    Helper Functions for mibrr
### Author:   Kyle M. Lang
### Created:  2014-DEC-09
### Modified: 2016-NOV-08

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



## Print startup message:
.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file = system.file("DESCRIPTION", package = pkgname),
                        fields = "Version")
    
    greet <-
        strwrap(
            paste0("Loading: ",
                   pkgname,
                   " ",
                   version,
                   ", Copyright (C) 2016 Kyle M. Lang. ",
                   pkgname,
                   " is distributed under Version 3 of the GNU Lesser General",
                   " Public License (LGPL-3); execute 'mibrrL()' for details. ",
                   pkgname,
                   " comes with ABSOLUTELY NO WARRANTY; execute 'mibrrW()' for",
                   " details. ",
                   pkgname,
                   " is beta software. Please report any bugs."),
            width = 81)
    
    for(i in greet) packageStartupMessage(i)
}



## Calculate the potential scale reduction factor (R-Hat)
calcRHat <- function(simsIn, nChains = 1)
{
    subChainLen <- floor(length(simsIn) / 2)
    nSubChains  <- nChains * 2

    if(length(simsIn) %% nSubChains == 0) {
        simsMat <- matrix(simsIn, ncol = nSubChains)
    } else {
        simsMat <- matrix(
            simsIn[1 : (length(simsIn) - (nSubChains - 1))],
            ncol = nSubChains
        )
    }

    wMean     <- colMeans(simsMat)
    grandMean <- mean(simsMat)

    bVar <- (subChainLen / (nSubChains - 1)) * sum((wMean - grandMean)^2)
    wVar <- mean(apply(simsMat, 2, var))
    tVar <- ((subChainLen - 1) / subChainLen) * wVar + (1 / subChainLen) * bVar
    
    sqrt(tVar / wVar)
}# END calcRHat()



## Sample imputations from the target's posterior gibbs sample:
sampleImps <- function(x, nDraws, nImps)
{
    sampledIndex <- sample(c(1 : nDraws), nImps)
    x$imps[sampledIndex, ]
}



## Fill the missing data with the sampled imputations:
fillMissing <- function(impNum,
                        targetVars,
                        targetMeans,
                        impSams,
                        data,
                        missList)
{
    ## Set number of previous generations to 3 to get back to the execution
    ## environment of mibrr():
    env <- parent.frame(3)
    
    for(j in targetVars)
        data[missList[[j]], j] <- impSams[[j]][impNum, ] + targetMeans[j]
    
    ## Restructure imputed data to match the raw data layout:
    if(env$control$preserveStructure)
        data <- data.frame(data, env$ignoredColumns)[ , env$rawNames]
    data
}



## Sample the imputations from the stationary posterior predictive distibution
## of the missing data
getImputedData <- function(gibbsState,
                           nImps,
                           data,
                           targetVars,
                           targetMeans,
                           missList)
{
    nDraws <- nrow(gibbsState[[1]]$imps)
    
    impSams <- lapply(X      = gibbsState,
                      FUN    = sampleImps,
                      nDraws = nDraws,
                      nImps  = nImps)

    lapply(X           = c(1 : nImps),
           FUN         = fillMissing,
           targetVars  = targetVars,
           targetMeans = targetMeans,
           impSams     = impSams,
           data        = data,
           missList    = missList)
}# END getImputedData()



## Use a variant of the method recommended by Park and Casella (2008) to get
## starting values for the MIBL penalty parameters
getLambdaStarts <- function(inData, nTargets, nPreds, nSamples = 25)
{
    if(require(mice)) {# Can we load mice()?
        
        ## Fill any missing data with rough guesses:
        micePreds <- quickpred(inData)
        miceOut <- mice(data            = inData,
                        m               = 1,
                        method          = "norm",
                        predictorMatrix = micePreds,
                        printFlag       = FALSE)
        
        impData <- as.matrix(complete(miceOut, 1))
        lambdaStart <- vector("numeric", nTargets)
        
        for(i in 1 : nTargets) {
            if(0.90 * nrow(impData) > (ncol(impData) - 1)) {# P << N
                tmpPredCount   <- ncol(impData)
                tmpOut         <- lm(impData[ , i] ~ impData[ , -i])
                lambdaStart[i] <-
                    tmpPredCount * sqrt(anova(tmpOut)["Residuals", "Mean Sq"]) /
                        sum(abs(tmpOut$coefficients[-1]))
            } else {
                ## If P ~ N or  P > N, subsample inData's columns
                ## and repeatedly apply the Park & Casella (2008) method.
                tmpLambda    <- vector("numeric", nSamples)
                tmpPredCount <- round(0.90 * nrow(impData))
                for(j in 1 : nSamples) {
                    predSelector <-
                        sample(c(1 : ncol(impData))[-i], size = tmpPredCount)
                    tmpDat       <- impData[ , predSelector]
                    tmpOut       <- lm(impData[ , i] ~ tmpDat)
                    tmpLambda[j] <- tmpPredCount *
                        sqrt(anova(tmpOut)["Residuals", "Mean Sq"]) /
                            sum(abs(tmpOut$coefficients[-1]))
                }# END for(j in 1 : nSamples)
                lambdaStart[i] <- mean(tmpLambda)
            }# END if( nrow(inData) > ncol(inData) )     
        }# END for(i in 1 : nTargets)
        
    } else {# mice() isn't available    
        
        warning("I cannot attach the mice package, so I cannot give lambda \
starting values via\nthe Park and Casella (2008) rule. I am falling back to the \
default values.")
        lambdaStarts <- rep(0.5, nTargets)  
        
    }# END if(require(mice))
    lambdaStarts
}# END getLambdaStarts()
    


## Compute R-Hat values for model parameters and throw a warning if any R-Hats
## exceed a given threshold
checkGibbsConv <- function(targetName,
                           gibbsStates,
                           critVal,
                           returnRHats = TRUE)
{
    gibbsState <- gibbsStates[[targetName]]
    
    ## Compute R-Hat values to check convergence:
    betaRHats <- apply(gibbsState$beta, 2, FUN = calcRHat, nChains = 1)
    tauRHats  <- apply(gibbsState$tau, 2, FUN = calcRHat, nChains = 1)
    sigmaRHat <- calcRHat(gibbsState$sigma, nChains = 1)
    
    ## Find nonconvergent Gibbs samples:
    badBetaCount <- sum(betaRHats > critVal)
    maxBetaRHat  <- max(betaRHats)
    badTauCount  <- sum(tauRHats > critVal)
    maxTauRHat   <- max(tauRHats)
    badSigmaFlag <- sigmaRHat > critVal

    ## Return warnings about possible failures of convergence:
    if(badBetaCount > 0) {
        warning(paste0("While imputing ",
                       targetName,
                       ", Beta's final Gibbs sample may not have converged.\n",
                       badBetaCount,
                       " R-Hats > ",
                       critVal,
                       " with maximum R-Hat = ",
                       round(maxBetaRHat, 4),
                       ".\nConsider increasing the size of the ",
                       "(retained) Gibbs samples."))
    }
    if(badTauCount > 0) {
        warning(paste0("While imputing ",
                       targetName,
                       ", Tau's final Gibbs sample may not have converged.\n",
                       badTauCount,
                       " R-Hats > ",
                       critVal,
                       " with maximum R-Hat = ",
                       round(maxTauRHat, 4),
                       ".\nConsider increasing the size of the ",
                       "(retained) Gibbs samples."))
    }
    if(badSigmaFlag) {
        warning(paste0("While imputing ",
                       targetName,
                       ", Sigma's final Gibbs sample ",
                       "may not have converged.\nR-Hat = ",
                       round(sigmaRHat, 4),
                       ".\nConsider increasing the size of the ",
                       "(retained) Gibbs samples."))
    }
    if(returnRHats) {
        list(beta = betaRHats, tau = tauRHats, sigma = sigmaRHat)
    } else {
        list()
    }
}# END checkGibbsConv()



## Initialize the gibbs sampled parameters with draws from the parameters'
## respective prior distributions
initializeParams <- function(data, nTargets, doBl, control)
{    
    nRows       <- nrow(data)
    nObsVec     <- colSums(!is.na(data))
    nPreds      <- ncol(data) - 1
    dataScales  <- apply(data, 2, FUN = sd, na.rm = TRUE)
    sigmaStarts <- dataScales[1 : nTargets]
    tauStarts   <- betaStarts <- matrix(NA, nPreds, nTargets)
    dvStarts    <- matrix(NA, nRows, nTargets)

    ## Populate the starting values for Lambda:
    if(!doBl) {
        options(warn = -1)# Suppress warnings about recycling elements
        lambda1Starts <- matrix(control$lambda1Starts, nTargets, 1)
        lambda2Starts <- matrix(control$lambda2Starts, nTargets, 1)
        options(warn = 0)
        lambdaMat <- cbind(lambda1Starts, lambda2Starts)
    } else {
        if(control$usePCStarts) {
            ## Must call this before the NA's are replaced with missCode:
            lambdaVec <- getLambdaStarts(inData   = data,
                                         nTargets = nTargets,
                                         nPreds   = nPreds)
        } else {
            options(warn = -1)
            lambdaVec <- as.vector(matrix(control$lambda1Starts, nTargets, 1))
            options(warn = 0)
        }
        lambdaMat <- cbind(lambdaVec, 0)
    }# END if(!doBl)
    
    ## Populate starting values for betas, taus, and sigma:
    for(j in 1 : nTargets) {
        if(!doBl) {
            lam1 <- lambdaMat[j, 1]
            lam2 <- lambdaMat[j, 2]
            
            tauPriorScale <- (8 * lam2 * sigmaStarts[j]) / lam1^2
            
            for(k in 1 : nPreds) {
                tauDraw <- 0.0
                while(tauDraw < 1.0)
                    tauDraw <- rgamma(1, shape = 0.5, scale = tauPriorScale)
                tauStarts[k, j] <- tauDraw
            }
            
            betaPriorCov <- diag(1 / (
                (lam2 / sigmaStarts[j]) *
                    (tauStarts[ , j] / (tauStarts[ , j] - 1.0))
            ))
            
            betaStarts[ , j] <- rmvnorm(1, rep(0, nPreds), betaPriorCov)
        } else {# We're doing MIBL
            lam              <- lambdaMat[j, 1]
            tauStarts[ , j]  <- rexp(nPreds, rate = (0.5 * lam^2))
            betaPriorCov     <- sigmaStarts[j] * diag(tauStarts[ , j])
            betaStarts[ , j] <- rmvnorm(1, rep(0, nPreds), betaPriorCov)
        }
    }# CLOSE for(j in 1 : nTargets)
       
    list(lambda = lambdaMat,
         beta   = betaStarts,
         tau    = tauStarts,
         sigma  = sigmaStarts)
}# END initializeParams()



## Compute the posterior means of the Gibbs samples
gibbsMeanFun <- function(gibbsState)
{
    outList <- list()
    outList$sigma <- mean(gibbsState$sigma)
    outList$tau   <- colMeans(gibbsState$tau)
    outList$beta  <- colMeans(gibbsState$beta)
    outList$dv    <- colMeans(gibbsState$imputedVar)
    
    outList
}



## Make sure that all control parameters are initialized
padControlList <- function()
{
    env <- parent.frame()
    ## Define the default control parameters:
    defaults = list(
        approxBurn        = env$sampleSizes[1],
        tuneBurn          = env$sampleSizes[2],
        postBurn          = env$sampleSizes[3],
        convThresh        = 1.1,
        lambda1Starts     = rep(0.5, env$nTargets),
        lambda2Starts     = rep(env$nPreds / 10, env$nTargets),
        usePCStarts       = FALSE,
        mcemEpsilon       = 1.0e-5,
        smoothingWindow   = 1,
        center            = TRUE,
        scale             = TRUE,
        adaptScales       = TRUE,
        simpleIntercept   = FALSE,
        twoPhaseOpt       = TRUE,
        minPredCor        = 0.3,
        miceIters         = 10,
        miceRidge         = 1e-4,
        miceMethod        = "pmm",
        fimlStarts        = FALSE,
        preserveStructure = TRUE
    )
    
    ## Pad the user-provided control list with default values:
    defaults[names(defaults) %in% names(env$control)] <- env$control
    defaults
}# END padControlList()



## Fill missing data with an appropriate integer code:
applyMissCode <- function() {
    env <- parent.frame()
    ## Construct an integer-valued missing data code that
    ## does not take legal data values and use it to flag NAs.
    if(is.null(env$missCode)) {
        if(max(abs(env$data), na.rm = TRUE) < 1.0) {
            env$missCode <- -9
        } else {
            codeMag <- floor(log10(max(abs(env$data), na.rm = TRUE))) + 2
            env$missCode <- -(10^codeMag - 1)
        }
    }
    env$data[is.na(env$data)] <- env$missCode
}



## Check the user inputs and isolate a set of target variables:
checkInputs <- function() {
    env <- parent.frame()
    
    ## Check for target variables. When no targets are given, all incomplete
    ## variables not listed in 'ignoreVars' are imputed.
    if(is.null(env$targetVars)) {
        if(env$doImp) {
            targetCandidates <-
                colnames(env$data)[!colnames(env$data) %in% env$ignoreVars]
            warning("You did not specify any target variables, so I will impute \
the missing data on\nevery variable in 'data' that is not listed in \
'ignoreVars'.\n")        
        } else {
            stop("Please specify a DV.")
        }
    } else {
        targetCandidates <- env$targetVars
    }
        
    ## Make sure 'data' contains missing data that we can find:
    if(env$doImp) {
        if(is.null(env$missCode)) {
            if(length(targetCandidates) > 1) {
                completeTargets <-
                    colMeans(is.na(env$data[ , targetCandidates])) == 0
            } else {
                completeTargets <-
                    mean(is.na(env$data[ , targetCandidates])) == 0
            }
            if(all(completeTargets)) {
                stop("Your target variables appear to be fully observed. Did \
you forget to provide a\nvalue for 'missCode'?\n")
            }
        } else {
            rMat <- env$data == env$missCode
            if(!any(rMat, na.rm = TRUE)) {
                stop(paste0("The value you provided for 'missCode' (i.e., ",
                            env$missCode,
                            ") does not appear anywhere in 'data'.\n",
                            "Are you sure that ",
                            env$missCode,
                            " encodes your missing data?\n")
                     )
            } else {
                env$data[rMat] <- NA
            }
        }
    }
    
    ## Select the final set of target variables:
    if(env$doImp) {
        env$targetVars <- targetCandidates[!completeTargets]
        if(any(completeTargets)) {
            warning(paste0("The potential target variables {",
                           paste(targetCandidates[completeTargets],
                                 collapse = ", "),
                           "} are fully observed.\n",
                           "These items will not be imputed.\n")
                    )
        }
    }
}# END checkInputs()



scaleDataWithFiml <- function(revert = FALSE) {
    env  <- parent.frame()
    nObs <- nrow(env$data)
    nVar <- ncol(env$data)

    if(!revert) {# Doing initial scaling
        ## Specify a lavaan model to estimate data's sufficient statistics:
        mod1 <- paste(
            paste0("F",
                   colnames(env$data),
                   " =~ 1*",
                   colnames(env$data),
                   "\n"),
            collapse = "")
        
        ## Estimate the sufficient statistics with FIML:
        out1 <- lavaan(model           = mod1,
                       data            = env$data,
                       int.ov.free     = FALSE,
                       int.lv.free     = TRUE,
                       auto.var        = TRUE,
                       auto.fix.single = TRUE,
                       missing         = "fiml")
        
        ## Store the item means and scales:
        env$dataMeans  <- as.vector(inspect(out1, "coef")$alpha)
        if(env$control$scale)
            env$dataScales <- sqrt(diag(inspect(out1, "coef")$psi))
        else
            env$dataScales <- rep(1, nVar)
        
        names(env$dataMeans) <- names(env$dataScales) <- colnames(env$data)
        
        ## Mean center data:
        if(env$control$center)
            env$data <- env$data - data.frame(
                matrix(env$dataMeans, nObs, nVar, byrow = TRUE)
            )
    } else {# Reverting the data to its original scaling
        env$data <- env$data + data.frame(
            matrix(env$dataMeans, nObs, nVar, byrow = TRUE)
        )
    }
        
}# END scaleData()



scaleData <- function(revert = FALSE) {
    env  <- parent.frame()
    nObs <- nrow(env$data)
    nVar <- ncol(env$data)

    if(!revert) {# Doing initial scaling
        if(env$control$scale)
            env$dataScales <- unlist(lapply(env$data, sd))
        else
            env$dataScales <- rep(1, nVar)
        
        ## Mean center data:
        if(env$control$center) {
            env$dataMeans <- colMeans(env$data)
            env$data   <- as.data.frame(
                scale(env$data, center = TRUE, scale = FALSE)
            )
        } else {
            env$dataMeans <- rep(0, nVar)
        }
        
        names(env$dataMeans) <- names(env$dataScales) <- colnames(env$data)
        
    } else {# Reverting the data to its original scaling
        env$data <- env$data + data.frame(
            matrix(env$dataMeans, nObs, nVar, byrow = TRUE)
        )
    }
}# END scaleData()



imputeCovs <- function() {
    env      <- parent.frame()
    covNames <- setdiff(colnames(env$data), env$targetVars)
    
    ## Construct a predictor matrix for mice() to use:
    predMat <- quickpred(env$data, mincor  = env$control$minPredCor)

    ## Construct a vector of elementary imputation methods:
    methVec           <- rep("", ncol(env$data))
    names(methVec)    <- colnames(env$data)
    methVec[covNames] <- env$control$miceMethod

    ## Singly impute the missing covariate values:
    miceOut <- mice(data            = env$data,
                    m               = 1,
                    maxit           = env$control$miceIters,
                    method          = methVec,
                    predictorMatrix = predMat,
                    printFlag       = FALSE,
                    ridge           = env$control$miceRidge)
    
    ## Replace missing covariate values with their imputations:
    env$data[ , covNames] <- complete(miceOut)[ , covNames]
}# END imputeCovs()



## Initially fill the missing values via single imputation:
simpleImpute <- function() {
    env  <- parent.frame()
    rFlags <- colMeans(is.na(env$data)) > 0
   
    ## Construct a predictor matrix for mice() to use:
    predMat <- quickpred(env$data, mincor = env$control$minPredCor)
    
    ## Construct a vector of elementary imputation methods:
    methVec       <- rep("", ncol(env$data))
    methVec[rFlags] <- env$control$miceMethod

    ## Singly impute the missing covariate values:
    miceOut <- mice(data            = env$data,
                    m               = 1,
                    maxit           = env$control$miceIters,
                    method          = methVec,
                    predictorMatrix = predMat,
                    printFlag       = FALSE,
                    ridge           = env$control$miceRidge)
    
    ## Replace missing values with their imputations:
    env$data <- complete(miceOut)
}# END imputeCovs()



nameOutput <- function() {
    env <- parent.frame()
    
    if(env$returnConvInfo) {
        names(env$rHatList) <- env$targetVars

        if(ncol(env$gibbsOut[[1]]$lambdaHistory) == 2)
            lamNames <- c("lambda1", "lambda2")
        else
            lamNames <- "lambda"
        
        for(v in env$targetVars)
            colnames(env$gibbsOut[[v]]$lambdaHistory) <- lamNames
    }
    
    if(env$returnParams)
        for(v in env$targetVars) {
            colnames(env$gibbsOut[[v]]$beta) <-
                c("intercept", setdiff(colnames(env$data), v))
            colnames(env$gibbsOut[[v]]$tau) <- setdiff(colnames(env$data), v)
        }
}



predictMibrr <- function(object,
                         newData,
                         targetVar,
                         nDraws = 0)
{
    if(nDraws == 0) {
        beta  <- matrix(colMeans(object$params[[targetVar]]$beta))
        sigma <- mean(object$params[[targetVar]]$sigma)
        
        out <- cbind(1, newData) %*% beta + rnorm(1, 0, sqrt(sigma))
    } else if(nDraws > 0) {
        index <- sample(c(1 : length(object$params[[targetVar]]$sigma)), nDraws)
        beta  <- object$params[[targetVar]]$beta[index, ]
        sigma <- object$params[[targetVar]]$sigma[index]
        
        out <- matrix(NA, nrow(newData), nDraws)
        for(j in 1 : nDraws)
            out[ , j] <- cbind(1, newData) %*% matrix(beta[j, ]) +
                rnorm(1, 0, sqrt(sigma[j]))
    } else {
        stop("nDraws must be non-negative.")
    }
    out
}



simulateData <- function(nObs,
                         nPreds,
                         r2,
                         collin,
                         beta,
                         means           = 0,
                         scales          = 1,
                         latentStructure = FALSE,
                         itemsPerFactor  = 1,
                         itemReliability = NULL)
{
    if(length(means) == 1) means <- rep(means, nPreds)
    
    w1 <- matrix(scales, nPreds, nPreds)
    w2 <- matrix(scales, nPreds, nPreds, byrow = TRUE)
    
    maxCov <- w1*w2
    
    sigma <- maxCov * collin
    diag(sigma) <- scales^2
    
    X <- cbind(1, rmvnorm(nObs, means, sigma))

    eta <- X %*% beta
    sigmaY <- (var(eta) / r2) - var(eta)
    y <- eta + rnorm(nObs, 0, sqrt(sigmaY))

    if(latentStructure) {
        nItems   <- nPreds * itemsPerFactor
        loadings <- matrix(0, nItems, nPreds)
        
        for(m in 1 : nPreds) {
            for(n in 1 : itemsPerFactor) {
                offset <- (m - 1) * itemsPerFactor
                loadings[n + offset, m] <- sqrt(itemReliability)
            }
        }

        theta <- diag(rep(1 - itemReliability, nItems))
    
        X <- X[ , -1] %*% t(loadings) + rmvnorm(nObs, rep(0, nItems), theta)
        
        outDat <- data.frame(y, X)
        colnames(outDat) <- c("y", paste0("x", c(1 : nItems)))
    } else {
        outDat <- data.frame(y, X[ , -1])
        colnames(outDat) <- c("y", paste0("x", c(1 : nPreds)))
    }
    outDat
}
