### Title:    Optimization Functions for MIBRR
### Author:   Kyle M. Lang
### Created:  2017-SEP-30
### Modified: 2017-NOV-28

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
        log(pgamma(l1^2 / (8 * sigmas * l2), 0.5, lower.tail = FALSE) *
            gamma(0.5)
            )
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
    (1 / (pgamma(tmp, 0.5, lower.tail = FALSE) * gamma(0.5))) *
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
