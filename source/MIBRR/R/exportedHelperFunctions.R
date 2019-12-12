### Title:    Exported Helper Functions for MIBRR
### Author:   Kyle M. Lang
### Created:  2014-12-09
### Modified: 2019-12-12

##--------------------- COPYRIGHT & LICENSING INFORMATION --------------------##
##  Copyright (C) 2019 Kyle M. Lang <k.m.lang@uvt.nl>                         ##
##                                                                            ##
##  This file is part of MIBRR.                                               ##
##                                                                            ##
##  This program is free software: you can redistribute it and/or modify it   ##
##  under the terms of the GNU General Public License as published by the     ##
##  Free Software Foundation, either version 3 of the License, or (at you     ##
##  option) any later version.                                                ##
##                                                                            ##
##  This program is distributed in the hope that it will be useful, but       ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of                ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General  ##
##  Public License for more details.                                          ##
##                                                                            ##
##  You should have received a copy of the GNU General Public License along   ##
##  with this program. If not, see <http://www.gnu.org/licenses/>.            ##
##----------------------------------------------------------------------------##


## Sample the imputations from the stationary posterior predictive distibution
## of the missing data
getImpData <- function(mibrrFit, nImps) {
    impList <- list()
    for(m in 1 : nImps) impList[[m]] <- mibrrFit$getImpDataset()
    impList
}# END getImputedData()

###--------------------------------------------------------------------------###

## Extract the parameter samples from a fitted MibrrFit object:
getParams <- function(mibrrFit, target, mix = TRUE) {
    ## Define an appropriate extractor function:
    fun <-
        switch(as.numeric(mix) + 1, mibrrFit$getSamples, mibrrFit$poolSamples)
    
    ## Populate the output object:
    out       <- list()
    out$beta  <- fun("beta", target)
    out$sigma <- fun("sigma", target)
    
    if(mibrrFit$penalty != 0) {# Used shrinkage priors?
        out$tau <- fun("tau", target)

        ## Don't mix MCEM chains:
        if(mibrrFit$doMcem)
            fun <- mibrrFit$getSamples
        
        out$lambda1 <- fun("lambda1", target)
        
        if(mibrrFit$penalty == 2)
            out$lambda2 <- fun("lambda2", target)
    }
    out
}

###--------------------------------------------------------------------------###

## Extract the posterior predictive samples from a fitted MibrrFit object:
getPostPredSams <- function(mibrrFit, target, mix = TRUE) {
    ## Define an appropriate extractor function:
    fun <-
        switch(as.numeric(mix) + 1, mibrrFit$getSamples, mibrrFit$poolSamples)
    
    ## Extract and return the posterior samples:
    fun("ppSams", target)
}

###--------------------------------------------------------------------------###

## Generate posterior predictions from a fitted BEN or BL model:
postPredict <- function(mibrrFit,
                        newData,
                        targetVars = NULL,
                        nDraws     = 0,
                        scale      = FALSE,
                        cenType    = "mode")
{
    if(!is.data.frame(newData)) stop("'newData' must be a data.frame")
    if(is.null(targetVars))     targetVars <- mibrrFit$targetVars
    
    if(scale) newData <- scale(newData)
    
    outList <- list()
    for(nm in targetVars) {
        pars     <- getParams(mibrrFit, nm)
        testData <- cbind(1, as.matrix(newData[ , colnames(pars$tau)]))
        
        if(nDraws > 0) {
            index <- sample(c(1 : length(pars$sigma)), nDraws)
            beta  <- pars$beta[index, ]
            sigma <- pars$sigma[index]
            
            out <- matrix(NA, nrow(testData), nDraws)
            for(j in 1 : nDraws)
                out[ , j] <- testData %*% matrix(beta[j, ]) +
                    rnorm(1, 0, sqrt(sigma[j]))
        }
        else {
            if(nDraws == 0) {
                beta  <- matrix(apply(pars$beta, 2, numMode))
                sigma <- numMode(pars$sigma)
            }
            else {# nDraws < 0
                cenTen <-
                    switch(cenType,
                           mode   = numMode,
                           median = median,
                           mean   = mean,
                           stop("'cenType' must be one of {'mode', 'median', 'mean'}")
                           )
                
                beta  <- matrix(apply(pars$beta, 2, cenTen))
                sigma <- cenTen(pars$sigma)
            }
            out <- testData %*% beta + rnorm(1, 0, sqrt(sigma))
        }
        outList[[nm]] <- out
    }
    outList
}

###--------------------------------------------------------------------------###

## Plot the posterior predictive density of targetVars vs. their observed
## densities:
ppCheck <- function(mibrrFit, targetVars = NULL, nSams = NULL) {
    
    if(!mibrrFit$control$savePpSams)
        stop("The object provided for 'mibrrFit' does not contain posterior predictive samples.")
    
    if(is.null(targetVars)) targetVars <- mibrrFit$targetVars

    ## Extract the posterior predictive samples:
    ppSams0        <- lapply(X   = targetVars,
                             FUN = function(x, y) getPostPredSams(y, x),
                             y   = mibrrFit)
    names(ppSams0) <- targetVars

    ## Define the index of rows to sample:
    index <- 1 : nrow(ppSams0[[1]])
    if(!is.null(nSams)) {
        check <- nSams > length(index)
        if(check) {
            warning(
                paste0("You have requested too many samples, so I will use all of the ",
                       length(index),
                       " samples stored in the provided 'mibrrFit' object.")
            )
            nSams <- length(index)
        }
        index <- sample(index, nSams)
    }

    ## Plot the posterior predictive densities:
    for(v in targetVars) {
        ppSams <- ppSams0[[v]][index, ]
        
        d0 <- density(mibrrFit$data[ , v], na.rm = TRUE)
        d1 <- apply(ppSams, 1, density)
        
        plot(NULL,
             ylim = range(d0$y, lapply(d1, "[[", x = "y")),
             xlim = range(d0$x, lapply(d1, "[[", x = "x")),
             main = paste0("Variable = ",
                           v,
                           "\nPP Densities (Red) vs. Obs. Density (Black)"),
             ylab = "Density",
             xlab = v)
        
        lapply(d1, lines, col = "red")
        lines(d0, col = "black", lwd = 2.5)
    }
}

###--------------------------------------------------------------------------###
    
## Access arbitrary fields in a 'MibrrFit' object:
getField <- function(mibrrFit, what) mibrrFit$field(what)
