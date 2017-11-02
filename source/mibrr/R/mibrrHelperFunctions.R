### Title:    Helper Functions for mibrr
### Author:   Kyle M. Lang
### Created:  2014-DEC-09
### Modified: 2017-NOV-02

##--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------##
##  Copyright (C) 2017 Kyle M. Lang <kyle.lang@ttu.edu>                        ##  
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
    version <- read.dcf(file   = system.file("DESCRIPTION", package = pkgname),
                        fields = "Version")
    
    greet <-
        strwrap(
            paste0("Loading: ",
                   pkgname,
                   " ",
                   version,
                   ", Copyright (C) ",
                   format(Sys.time(), "%Y"),
                   " Kyle M. Lang. ",
                   pkgname,
                   " is distributed under Version 3 of the GNU Lesser General Public License (LGPL-3); execute 'mibrrL()' for details. ",
                   pkgname,
                   " comes with ABSOLUTELY NO WARRANTY; execute 'mibrrW()' for details. ",
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
        miceOut   <- mice(data            = inData,
                          m               = 1,
                          method          = "norm",
                          predictorMatrix = micePreds,
                          printFlag       = FALSE)
        
        impData     <- as.matrix(complete(miceOut, 1))
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
        
        warning("I cannot attach the mice package, so I cannot give lambda starting values via\nthe Park and Casella (2008) rule. I am falling back to the default values.")
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

    ## NOTE: We don't need to start the intercept. It's initial value will be
    ##       sampled in the first iteration of the Gibbs sampler.
    
    ## Populate the starting values for Lambda:
    if(!doBl) {
        options(warn = -1)# Suppress warnings about recycling elements
        lambda1Starts <- matrix(control$lambda1Starts, nTargets, 1)
        lambda2Starts <- matrix(control$lambda2Starts, nTargets, 1)
        options(warn = 0)
        lambdaMat <- cbind(lambda1Starts, lambda2Starts)
    } else {
        if(control$usePcStarts) {
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
        } else {# We're doing BL
            lam <- lambdaMat[j, 1]
            
            tauStarts[ , j] <- rexp(nPreds, rate = (0.5 * lam^2))

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
        convThresh        = 1.1,
        lambda1Starts     = rep(0.5, env$nTargets),
        lambda2Starts     = rep(env$nPreds / 10, env$nTargets),
        usePcStarts       = FALSE,
        smoothingWindow   = 1,
        center            = TRUE,
        scale             = TRUE,
        adaptScales       = TRUE,
        simpleIntercept   = FALSE,
        minPredCor        = 0.3,
        miceIters         = 10,
        miceRidge         = 1e-4,
        miceMethod        = "pmm",
        fimlStarts        = FALSE,
        preserveStructure = TRUE,
        optTraceLevel     = 0,
        optCheckKkt       = TRUE,
        optMethod         = "L-BFGS-B",
        optBoundLambda    = TRUE #,
        #optReturnACov     = FALSE
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
            warning("You did not specify any target variables, so I will impute the missing data on\nevery variable in 'data' that is not listed in 'ignoreVars'.\n")        
        } else {
            stop("Please specify a DV.")
        }
    } else {
        targetCandidates <- env$targetVars
    }
        
    ## Make sure 'data' contains missing data that we can find:
                                        #if(env$doImp) {
    if(is.null(env$missCode)) {
        rMat <- is.na(env$data)
    } else {
        rMat <- env$data == env$missCode
        
        if(!any(rMat, na.rm = TRUE))
            stop(paste0("The value you provided for 'missCode' (i.e., ",
                        env$missCode,
                        ") does not appear anywhere in 'data'.\nAre you sure that ",
                        env$missCode,
                        " encodes your missing data?\n")
                 )
    }
        
    if(length(targetCandidates) > 1) 
        completeTargets <- colMeans(rMat[ , targetCandidates]) == 0
    else 
        completeTargets <- mean(rMat[ , targetCandidates]) == 0
    
    if(env$doImp & all(completeTargets)) 
        stop("Your target variables appear to be fully observed. Did you forget to provide a\nvalue for 'missCode'?\n")
    
    ## Select the final set of target variables:
    if(env$doImp) {
        env$targetVars <- targetCandidates[!completeTargets]
        if(any(completeTargets))
            warning(
                paste0("The potential target variables {",
                       paste(targetCandidates[completeTargets], collapse = ", "),
                       "} are fully observed.\nThese items will not be imputed.\n")
            )
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
        
}# END scaleDataWithFiml()



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
            env$data      <- as.data.frame(
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

    ## Don't try to impute a fully observed covariate matrix:
    check <- all(!is.na(env$data[ , covNames]))
    if(check) return()
    
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
    env$data[ , covNames] <- complete(miceOut, 1)[ , covNames]
}# END imputeCovs()



## Initially fill the missing values via single imputation:
simpleImpute <- function() {
    env    <- parent.frame()
    rFlags <- colSums(is.na(env$data)) > 0
   
    ## Construct a predictor matrix for mice() to use:
    predMat <- quickpred(env$data, mincor = env$control$minPredCor)
    
    ## Construct a vector of elementary imputation methods:
    methVec         <- rep("", ncol(env$data))
    methVec[rFlags] <- env$control$miceMethod

    ## Singly impute the missing values:
    miceOut <- mice(data            = env$data,
                    m               = 1,
                    maxit           = env$control$miceIters,
                    method          = methVec,
                    predictorMatrix = predMat,
                    printFlag       = FALSE,
                    ridge           = env$control$miceRidge)
    
    ## Replace missing values with their imputations:
    env$data <- complete(miceOut, 1)
}# END imputeCovs()



nameOutput <- function() {
    env <- parent.frame()
    
    if(env$returnConvInfo) names(env$rHatList) <- env$targetVars
    
    if(env$returnParams)
        for(v in env$targetVars) {
            colnames(env$gibbsOut[[v]]$beta) <-
                c("intercept", setdiff(colnames(env$data), v))
            colnames(env$gibbsOut[[v]]$tau) <- setdiff(colnames(env$data), v)
        }
}



predictMibrr <- function(object,
                         newData,
                         targetVar = NULL,
                         nDraws    = 0)
{
    if(!is.null(targetVar)) object$params <- object$params[targetVar]
    
    outList <- list()
    outInd  <- 0
    for(obj in object$params) {
        outInd <- outInd + 1
        if(nDraws == 0) {
            beta  <- matrix(colMeans(obj$beta))
            sigma <- mean(obj$sigma)
            
            out <- cbind(1, newData) %*% beta + rnorm(1, 0, sqrt(sigma))
        } else if(nDraws > 0) {
            index <- sample(c(1 : length(obj$sigma)), nDraws)
            beta  <- obj$beta[index, ]
            sigma <- obj$sigma[index]
            
            out <- matrix(NA, nrow(newData), nDraws)
            for(j in 1 : nDraws)
                out[ , j] <- cbind(1, newData) %*% matrix(beta[j, ]) +
                    rnorm(1, 0, sqrt(sigma[j]))
        } else {
            stop("nDraws must be non-negative.")
        }
        outList[[outInd]] <- out
    }
    names(outList) <- names(object$params)
    outList
}


##### OPTIMIZATION FUNCTIONS #####


### The conditional loglikelihood function of Lambda for use during the
### empirical bayes updating.
eNetLL <- function(lambdaVec, gibbsState) {
    l1 <- lambdaVec[1]
    l2 <- lambdaVec[2]
    
    taus   <- gibbsState$tau
    sigmas <- gibbsState$sigma
    betas  <- gibbsState$beta
    
    p <- ncol(taus)
    
    e1 <- mean(
        log(pgamma(l1^2 / (8 * sigmas * l2), 0.5, lower = FALSE) * gamma(0.5))
    )
    e2 <- mean(rowSums((taus / (taus - 1)) * betas[ , -1]^2) / sigmas)
    e3 <- mean(rowSums(taus) / sigmas)
    
    p * log(l1) - p * e1 - (l2 / 2) * e2 - (l1^2 / (8 * l2)) * e3 # LL
}# END eNetLL()



### The gradient function for the conditional LL of Lambda:
eNetGrad <- function(lambdaVec, gibbsState)
{
    l1 <- lambdaVec[1]
    l2 <- lambdaVec[2]

    taus   <- gibbsState$tau
    sigmas <- gibbsState$sigma
    betas  <- gibbsState$beta
    
    p   <- ncol(taus)
    tmp <- l1^2 / (8 * sigmas * l2)
    
    e1 <- mean(
    (1 / (pgamma(tmp, 0.5, lower = FALSE) * gamma(0.5))) *
    (1 / (sqrt(tmp) * exp(tmp))) * (1 / sigmas)
    )
    e2 <- mean(rowSums((taus / (taus - 1)) * betas[ , -1]^2) / sigmas)
    e3 <- mean(rowSums(taus) / sigmas)
    
    w1 <- l1 / (4 * l2)
    w2 <- l1^2 / (8 * l2^2)
    
    c((p / l1) + (p * w1 * e1) - (w1 * e3),  # dLL / dl1
    (-p * w2 * e1) - (0.5 * e2) + (w2 * e3)) # dLL / dl2
}# END eNetGrad()



## Wrapper to allow optimx to run within lapply():
optWrap <- function(targetIndex,
                    lambdaMat,
                    optFun,
                    optGrad,
                    optMethod,
                    optLower,
                    optHessian,
                    optControl,
                    myGibbs)
{
    optOut <- optimx(par        = lambdaMat[targetIndex, ],
                     fn         = optFun,
                     gr         = optGrad,
                     method     = optMethod,
                     lower      = optLower,
                     hessian    = optHessian,
                     control    = optControl,
                     gibbsState = myGibbs[[targetIndex]])
    
    if(length(optMethod) > 1) optOut <- optOut[nrow(optOut), ]
    
    tmpList <- list()
    if(optHessian) {
        hessMat      <- attr(optOut, "details")[ , "nhatend"][[1]]
        tmpList$vcov <- solve(-hessMat)
    }
    
    if(optControl$kkt) tmpList$kktFlags <- c(optOut$kkt1, optOut$kkt2)
    
    tmpList$lambda <- c(optOut[[1]], optOut[[2]])
    tmpList
}# END optWrap()



## Optimize the penalty parameter for the BL using the rule given in Park &
## Casella (2008):
updateBlLambda <- function(gibbsState) {
    taus <- gibbsState$tau
    p    <- ncol(taus)
    
    sqrt((2 * p) / sum(colMeans(taus)))
}



## Optimize the BEN or BL penalty parameters:
optimizeLambda <- function(lambdaMat,
                           gibbsState,
                           doBl       = FALSE,
                           returnACov = FALSE,
                           controlParms)
{    
    ## Use simple update rule and return early when doing BL:
    if(doBl) return(lapply(gibbsState, updateBlLambda))
    
    optMethod <- controlParms$method
    useSeqOpt <- length(optMethod) > 1
    
    if(controlParms$boundLambda) {
        lowBounds <- c(1e-5, 1e-5)
        optMethod <- "L-BFGS-B"
    } else {
        lowBounds <- -Inf
    }
    
    options(warn = ifelse(controlParms$showWarns, 0, -1))

    if(.Platform$OS.type == "unix") nullFile <- "/dev/null"
    else                            nullFile <- "nul"
    
    sink(nullFile) # Don't show optimx output
    
    optList <- lapply(c(1 : nrow(lambdaMat)),
                      FUN        = optWrap,
                      lambdaMat  = lambdaMat,
                      optFun     = eNetLL,
                      optGrad    = eNetGrad,
                      optMethod  = optMethod,
                      optLower   = lowBounds,
                      optHessian = returnACov,
                      optControl =
                          list(trace     = controlParms$traceLevel,
                               maximize  = TRUE,
                               kkt       = controlParms$checkKkt,
                               follow.on = useSeqOpt),
                      myGibbs    = gibbsState)
    
    sink()
    options(warn = 0)
    
    optList
}# END optimizeLambda()
