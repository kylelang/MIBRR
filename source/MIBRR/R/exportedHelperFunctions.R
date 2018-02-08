### Title:    Exported Helper Functions for MIBRR
### Author:   Kyle M. Lang
### Created:  2014-DEC-09
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



## Sample the imputations from the stationary posterior predictive distibution
## of the missing data
complete <- function(mibrrFit, nImps) {
    impList <- list()
    for(m in 1 : nImps) impList[[m]] <- mibrrFit$getImpDataset()
    impList
}# END getImputedData()



predictMibrr <- function(object,
                         newData,
                         targetVar = NULL,
                         nDraws    = 0)
{
    if(!is.data.frame(newData)) stop("'newData' must be a data.frame")
    if(!is.null(targetVar))     object$params <- object$params[targetVar]
    
    outList <- list()
    for(nm in names(object$params)) {
        obj      <- object$params[[nm]]
        testData <- cbind(1, as.matrix(newData[ , colnames(obj$tau)]))
        
        if(nDraws == 0) {
            beta  <- matrix(colMeans(obj$beta))
            sigma <- mean(obj$sigma)
            
            out <- testData %*% beta + rnorm(1, 0, sqrt(sigma))
        } else if(nDraws > 0) {
            index <- sample(c(1 : length(obj$sigma)), nDraws)
            beta  <- obj$beta[index, ]
            sigma <- obj$sigma[index]
            
            out <- matrix(NA, nrow(testData), nDraws)
            for(j in 1 : nDraws)
                out[ , j] <- testData %*% matrix(beta[j, ]) +
                    rnorm(1, 0, sqrt(sigma[j]))
        } else {
            stop("nDraws must be non-negative.")
        }
        outList[[nm]] <- out
    }
    outList
}
