### Title:    Helper Functions for mibrr
### Author:   Kyle M. Lang
### Created:  2014-DEC-09
### Modified: 2018-FEB-13

##--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------##
##  Copyright (C) 2018 Kyle M. Lang <k.m.lang@uvt.nl>                          ##  
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
                   " is distributed under Version 3 of the GNU General Public License (GPL-3); execute 'mibrrL()' for details. ",
                   pkgname,
                   " comes with ABSOLUTELY NO WARRANTY; execute 'mibrrW()' for details. ",
                   pkgname,
                   " is beta software. Please report any bugs."),
            width = 81)
    
    for(i in greet) packageStartupMessage(i)
}



## Use a variant of the method recommended by Park and Casella (2008) to get
## starting values for the MIBL penalty parameters
getLambdaStarts <- function(object, nSamples = 25)
{
    data     <- object$getData()
    nTargets <- object$nTargets
    
    ## Fill any missing data with rough guesses:
    micePreds <- quickpred(data)
    miceOut   <- mice(data            = data,
                      m               = 1,
                      method          = "norm",
                      predictorMatrix = micePreds,
                      printFlag       = FALSE)
    
    impData      <- as.matrix(complete(miceOut, 1))
    lambdaStarts <- vector("numeric", nTargets)
    
    for(i in 1 : nTargets) {
        if(0.90 * nrow(impData) > (ncol(impData) - 1)) {# P << N
            tmpPredCount    <- ncol(impData)
            tmpOut          <- lm(impData[ , i] ~ impData[ , -i])
            lambdaStarts[i] <-
                tmpPredCount * sqrt(anova(tmpOut)["Residuals", "Mean Sq"]) /
                sum(abs(tmpOut$coefficients[-1]))
        } else {
            ## If P ~ N or  P > N, subsample data's columns
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
            lambdaStarts[i] <- mean(tmpLambda)
        }# END if( nrow(data) > ncol(data) )     
    }# END for(i in 1 : nTargets)
    
    lambdaStarts
}# END getLambdaStarts()
    


compStatsWithFiml <- function(object, revert = FALSE) {
    nObs    <- object$nObs()
    nVar    <- object$nVar()
    dNames  <- object$dataNames()
    data    <- object$getData()
    control <- object$getControl()
    
    ## Specify a lavaan model to estimate data's sufficient statistics:
    mod1 <- paste(paste0("F", dNames, " =~ 1*", dNames, "\n"), collapse = "")
    
    ## Estimate the sufficient statistics with FIML:
    out1 <- lavaan(model           = mod1,
                   data            = data,
                   int.ov.free     = FALSE,
                   int.lv.free     = TRUE,
                   auto.var        = TRUE,
                   auto.fix.single = TRUE,
                   missing         = "fiml")
    
    ## Store the item means:
    tmp        <- as.vector(inspect(out1, "coef")$alpha)
    names(tmp) <- dNames
    object$setDataMeans(tmp)
    
    ## Store the item scales:
    if(control$scale) {
        tmp        <- sqrt(diag(inspect(out1, "coef")$psi))
        names(tmp) <- dNames
        object$setDataScales(tmp)
    }
    else {
        tmp        <- rep(1.0, nVar)
        names(tmp) <- dNames
        object$setDataScales(tmp)
    }
}# END scaleDataWithFiml()



## Initially fill the missing values via single imputation:
simpleImpute <- function(object, covsOnly = FALSE) {
    cn      <- object$dataNames()
    data    <- object$data
    control <- object$getControl()
    
    rFlags <- (object$countMissing() > 0)[cn]
    
    if(covsOnly) {
        impTargets <- setdiff(cn, object$targets())
        rFlags <- rFlags & cn %in% impTargets 
    }
    else {
        impTargets <- cn
    }
    
    ## Don't try to impute fully observed targets:
    if(!any(rFlags)) return()
    
    ## Construct a predictor matrix for mice() to use:
    predMat <- quickpred(data, mincor = control$minPredCor)
    
    ## Construct a vector of elementary imputation methods:
    methVec         <- rep("", ncol(data))
    methVec[rFlags] <- control$miceMethod
    
    ## Singly impute the missing values:
    miceOut <- mice(data            = data,
                    m               = 1,
                    maxit           = control$miceIters,
                    method          = methVec,
                    predictorMatrix = predMat,
                    printFlag       = FALSE,
                    ridge           = control$miceRidge)
    
    ## Replace missing values with their imputations:
    object$setData(complete(miceOut, 1)[ , rFlags])
}# END simpleImpute()


## Print 'x' only if 'verbose = TRUE':
vcat <- function(x) if(parent.frame()$mibrrFit$verbose) cat(x)


## Estimate the mode of a continous vector:
estMode <- function(x) {
    dens <- density(x, na.rm = TRUE)
    dens$x[which.max(dens$y)]
}
