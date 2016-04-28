### Title:    Helper Functions for mibrr
### Author:   Kyle M. Lang
### Created:  2014-DEC-09
### Modified: 2015-FEB-23

###################### COPYRIGHT & LICENSING INFORMATION ########################
###    Copyright (C) 2015 Kyle M. Lang <kylelang@ku.edu>                      ###  
###                                                                           ###
###    This program is free software: you can redistribute it and/or modify   ###
###    it under the terms of the GNU General Public License as published by   ###
###    the Free Software Foundation, either version 3 of the License, or      ###
###    (at your option) any later version.                                    ###
###                                                                           ###
###    This program is distributed in the hope that it will be useful,        ###
###    but WITHOUT ANY WARRANTY; without even the implied warranty of         ###
###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          ###
###    GNU General Public License for more details.                           ###
###                                                                           ###
###    You should have received a copy of the GNU General Public License      ###
###    along with this program.  If not, see <http://www.gnu.org/licenses/>.  ###
#################################################################################


## Calculate the potential scale reduction factor (R-Hat)
calcRHat <- function(simsIn, nChains = 1)
{
    subChainLen <- floor(length(simsIn) / 2)
    nSubChains <- nChains * 2

    if(length(simsIn) %% nSubChains == 0) {
        simsMat <- matrix(simsIn, ncol = nSubChains)
    } else {
        simsMat <-
            matrix(
                simsIn[1 : ( length(simsIn) - (nSubChains - 1) )],
                ncol = nSubChains
            )
    }

    wMean <- colMeans(simsMat)
    grandMean <- mean(simsMat)

    bVar <- ( subChainLen / (nSubChains - 1) ) *
        sum( (wMean - grandMean)^2 )
    wVar <- mean(apply(simsMat, 2, var))
    tVar <- ((subChainLen - 1) / subChainLen) *
        wVar + (1 / subChainLen) * bVar

    rHat <- sqrt(tVar / wVar)

    rHat
}# END calcRHat()


## Sample the imputations from the stationary posterior
## predictive distibution of the missing data
getImputedData <- function(gibbsState, nImps, rawData,
                           targetVars, targetMeans)
{
    nObs <- nrow(rawData)
    frozenLabels <- colnames(rawData)
    nTargets <- length(targetVars)
    nDraws <- nrow(gibbsState[[1]]$imps)
    predData <-
        rawData[ , !colnames(rawData) %in% targetVars]
    
    impSams <-
        lapply(gibbsState,
               FUN = function(x, nDraws, nImps){
                   sampledIndex <- sample(c(1 : nDraws), nImps)
                   x$imps[sampledIndex, ]
               },
               nDraws = nDraws,
               nImps = nImps)
    
    impDatList <-
        lapply(c(1 : nImps),
               FUN = function(x, targetVars, impSams, rawData, inLabs)
                   {
                       for(j in 1 : length(targetVars)) {
                           myDV <- targetVars[j]
                           naFlag <- is.na(rawData[ , myDV])
                           rawData[naFlag, myDV] <-
                               impSams[[j]][x, naFlag] + targetMeans[x]
                       }
                       colnames(rawData) <- inLabs
                       rawData
                   },
               targetVars = targetVars,
               impSams = impSams,
               rawData = rawData,
               inLabs = frozenLabels)
    
    impDatList
}# END getImputedData()


## Use a variant of the method recommended by Park and Casella (2008)
## to get starting values for the MIBL penalty parameters
getLambdaStarts <- function(inData, nTargets, nPreds, nSamples = 25)
{
    if(require(mice)) {# Can we load mice()?
        
        ## Fill any missing data with rough guesses:
        micePreds <- quickpred(inData)
        miceOut <- mice(inData,
                        m = 1,
                        method = "norm",
                        predictorMatrix = micePreds,
                        printFlag = FALSE)
        
        impData <- as.matrix(complete(miceOut, 1))
        lambdaStart <- vector("numeric", nTargets)
        
        for(i in 1 : nTargets) {
            if( 0.90 * nrow(impData) > (ncol(impData) - 1) ) {# P << N
                tmpPredCount <- ncol(impData)
                tmpOut <- lm(impData[ , i] ~ impData[ , -i])
                lambdaStart[i] <- tmpPredCount *
                    sqrt(anova(tmpOut)["Residuals", "Mean Sq"]) /
                        sum(abs(tmpOut$coefficients[-1]))
            } else {
                ## If P ~ N or  P > N, subsample inData's columns
                ## and repeatedly apply the Park & Casella (2008) method.
                tmpLambda <- vector("numeric", nSamples)
                tmpPredCount <- round( 0.90 * nrow(impData) )
                for(j in 1 : nSamples) {
                    predSelector <- sample(c( 1 : ncol(impData) )[-i],
                                           size = tmpPredCount)
                    tmpDat <- impData[ , predSelector]
                    tmpOut <- lm(impData[ , i] ~ tmpDat)
                    tmpLambda[j] <- tmpPredCount *
                        sqrt(anova(tmpOut)["Residuals", "Mean Sq"]) /
                            sum(abs(tmpOut$coefficients[-1]))
                }# END for(j in 1 : nSamples)
                
                lambdaStart[i] <- mean(tmpLambda)
                
            }# END if( nrow(inData) > ncol(inData) )     
        }# END for(i in 1 : nTargets)
        
    } else {# mice() isn't available    

        warning(paste0("I can't attach the mice package, so I can't give ",
                       "lambda starting values via the Park and Casella ",
                       "(2008) rule. I'm falling back to the default ",
                       "starting values."))
        lambdaStarts <- rep((nPreds / 10), nTargets)  

    }# END if(require(mice))
    
    lambdaStarts
}# END getLambdaStarts()
    
    
## Compute R-Hat values for model parameters and
## throw a warning if any R-Hats exceed a given threshold
checkGibbsConv <- function(targetIndex,
                           gibbsStates,
                           returnRHats = TRUE,
                           targetNames,
                           critVal)
{
    gibbsState <- gibbsStates[[targetIndex]]
    targetName <- targetNames[targetIndex]
    
    ## Compute R-Hat values to check convergence:
    betaRHats <-
        apply(gibbsState$beta, 2,
              FUN = calcRHat, nChains = 1)
    tauRHats <-
        apply(gibbsState$tau, 2,
              FUN = calcRHat, nChains = 1)
    sigmaRHat <-
        calcRHat(gibbsState$sigma, nChains = 1)

    badBetaCount <- sum(betaRHats > critVal)
    maxBetaRHat <- max(betaRHats)
    badTauCount <- sum(tauRHats > critVal)
    maxTauRHat <- max(tauRHats)
    badSigmaFlag <- sigmaRHat > critVal

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
        list(beta = betaRHats,
             tau = tauRHats,
             sigma = sigmaRHat)
    } else {
        list()
    }

}# END checkGibbsConv()


## Initialize the gibbs sampled parameters with draws
## from the parameters' respective prior distributions
initializeParams <- function(rawData,
                             inLambdas,
                             nTargets,
                             doENET = TRUE)
{    
    nRows <- nrow(rawData)
    nObsVec <- colSums(!is.na(rawData))
    nPreds <- ncol(rawData) - 1
    dataScales <- apply(rawData, 2, FUN = sd, na.rm = TRUE)

    sigmaStarts <- dataScales[1 : nTargets]
    intStarts <- vector("numeric", nTargets)
    tauStarts <- slopeStarts <-
        matrix(NA, nPreds, nTargets)
    dvStarts <- matrix(NA, nRows, nTargets)

    for(j in 1 : nTargets) {

        if(doENET) {
            lam1 <- inLambdas[j, 1]
            lam2 <- inLambdas[j, 2]
   
            tauPriorScale <-
                (8 * lam2 * sigmaStarts[j]) / lam1^2
            
            for(k in 1 : nPreds) {
                tauDraw <- 0.0
                while(tauDraw < 1.0) {
                    tauDraw <-
                        rgamma(1, shape = 0.5,
                               scale = tauPriorScale)
                }
                tauStarts[k, j] <- tauDraw
            }
            
            slopePriorCov <- diag(
                1 / (
                    (lam2 / sigmaStarts[j]) *
                        ( tauStarts[ , j] /
                             (tauStarts[ , j] - 1.0) )
                )
            )
            
            slopeStarts[ , j] <-
                rmvnorm(1, rep(0, nPreds), slopePriorCov)
        } else {
            lam <- inLambdas[j]

            tauStarts[ , j] <-
                rexp(nPreds, rate = (0.5 * lam^2))
            
            slopePriorCov <- 
                sigmaStarts[j] * diag(tauStarts[ , j])
            
            slopeStarts[ , j] <- rmvnorm(1,
                                         rep(0, nPreds),
                                         slopePriorCov)
        }
        
        ## Intercept starts are drawn from an informative,
        ## flat prior with bounds = +/- SEM
        intBound <-
            dataScales[j] / sqrt(nObsVec[j])

        intStarts[j] <- runif(1, -intBound, intBound)

    }# CLOSE for(j in 1 : nTargets)

    betaStarts <- rbind(intStarts, slopeStarts)

    list(beta = betaStarts,
         tau = tauStarts,
         sigma = sigmaStarts)

}# END initializeParams()


## Compute the posterior means of the Gibbs samples
gibbsMeanFun <- function(gibbsState)
{
    outList <- list()
    outList$sigma <- mean(gibbsState$sigma)
    outList$tau <- colMeans(gibbsState$tau)
    outList$beta <- colMeans(gibbsState$beta)
    outList$dv <- colMeans(gibbsState$imputedVar)
    
    outList
}


## Make sure that all control parameters are initialized
padControlList <- function(x, inList, defaultList)
{
    if(is.null(inList[[x]])) {
        defaultList[[x]]
    } else {
        inList[[x]]
    }
}



##### EXTRA STUFF ######

### The following functions are depricated, because the
### optimization of the penalty parameters is now done
### within the C++ source via nlopt.

## The conditional loglikelihood function of Lambda
## for use during the empirical Bayes updating.
eNetCondLL <- function(lambdaVec, gibbsState)
{
    lambda1 <- lambdaVec[1]
    lambda2 <- lambdaVec[2]
    taus <- gibbsState$tau
    sigmas <- gibbsState$sigma
    betas <- gibbsState$beta
    nPreds <- ncol(taus)

    weight1 <- nPreds * log(lambda1)
    weight2 <- lambda2 / (2 * sigmas)
    weight3 <- lambda1^2 / (8 * sigmas * lambda2)

    term1 <-
        log( pgamma(weight3, 0.5, lower = FALSE) * gamma(0.5) )
    
    term2 <-
        rowSums( ( taus / (taus - 1) ) * betas[ , -1]^2 )

    llVec <-
        weight1 - nPreds * term1 -
            weight2 * term2 -
                weight3 * rowSums(taus)
    
    condLL <- mean(llVec)

    condLL
}# END eNetCondLL()


## The gradient function for the conditional LL of Lambda
eNetGrad <- function(lambdaVec, gibbsState)
{
    lambda1 <- lambdaVec[1]
    lambda2 <- lambdaVec[2]
    taus <- gibbsState$tau
    sigmas <- gibbsState$sigma
    betas <- gibbsState$beta
    nPreds <- ncol(taus)

    weight11 <- nPreds / lambda1
    weight12 <- (nPreds * lambda1) / (4.0 * lambda2)
    weight13 <- lambda1 / (4.0 * sigmas * lambda2)
    weight21 <- (nPreds * lambda1^2) / (8.0 * lambda2^2)
    weight22 <- 1 / (2 * sigmas)
    weight23 <- lambda1^2 / (8.0 * sigmas * lambda2^2)

    tmpVal <- lambda1^2 / (8.0 * sigmas * lambda2)
    
    term1 <-
        ( 1 / ( pgamma(tmpVal, 0.5, lower = FALSE) * gamma(0.5) ) ) *
            ( 1 / sqrt(tmpVal) ) * exp(-tmpVal) * (1 / sigmas)
    
    term2 <-
        rowSums( ( taus / (taus - 1) ) * betas[ , -1]^2 )

    lambda1Grad <-
        mean( weight11 + (weight12 * term1) -
                 (weight13 * rowSums(taus))
             )
    lambda2Grad <-
        mean( -(weight21 * term1) - (weight22 * term2) +
                 (weight23 * rowSums(taus))
             )

    c(lambda1Grad, lambda2Grad)
}# END eNetGrad()


## Function to optimize the penalty parameters via
## numerical maximization of eNetCondLL()
optimizeLambda <- function(lambdaMat,
                           gibbsState,
                           printFlag = TRUE,
                           returnVCOV = FALSE,
                           controlParms)
{
    optMethod <- controlParms$method
    useSeqOpt <- controlParms$useSeqOpt

    if(controlParms$boundLambda) {
        lowBounds <- c(0, 0)
        optMethod <- "L-BFGS-B"
    } else {
        lowBounds <- -Inf
    }

    useSeqOpt <- ifelse(length(optMethod) > 1, useSeqOpt, FALSE)

    options( warn = ifelse(controlParms$showWarns, 0, -1) )

    if(!printFlag) sink("/dev/null")

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
            optOut <- optimx(par = lambdaMat[targetIndex, ],
                             fn = optFun,
                             gr = optGrad,
                             method = optMethod,
                             lower = optLower,
                             hessian = optHessian,
                             control = optControl,
                             gibbsState = myGibbs[[targetIndex]])

            if(length(optMethod) > 1) {
                optOut <- optOut[nrow(optOut), ]
            }

            tmpList <- list()
            if(optHessian) {
                hessMat <- attr(optOut, "details")[ , "nhatend"][[1]]
                tmpList$vcov <- solve(-hessMat)
            }
            
            if(optControl$kkt) {
                tmpList$kktFlags <- c(optOut$kkt1, optOut$kkt2)
            }
            
            tmpList$lambda <- c(optOut[[1]], optOut[[2]])

            tmpList
        }

    optList <- lapply(c(1 : nrow(lambdaMat)),
                      FUN = optWrap,
                      lambdaMat = lambdaMat,
                      optFun = eNetCondLL,
                      optGrad = eNetGrad,
                      optMethod = optMethod,
                      optLower = lowBounds,
                      optHessian = returnVCOV,
                      optControl =
                          list(trace = controlParms$traceLevel,
                               maximize = TRUE,
                               kkt = controlParms$checkKKT,
                               follow.on = useSeqOpt),
                      myGibbs = gibbsState)
    
    if(!printFlag) sink()
    options(warn = 0)
    
    optList
}# END optimizeLambda()
