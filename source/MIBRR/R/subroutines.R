### Title:    Subroutines for the MIBRR Package
### Author:   Kyle M. Lang
### Created:  2017-NOV-28
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


init <- function(doBl,
                 doImp,
                 doMcem,
                 data,
                 targetVars,
                 ignoreVars,
                 iterations,
                 sampleSizes,
                 lam1PriorPar,
                 lam2PriorPar,
                 missCode,
                 verbose,
                 seed,
                 control)
{
    if(!is.list(sampleSizes)) sampleSizes <- list(sampleSizes)
    
    ## Initialize a new MibrrFit object:
    mibrrFit <- MibrrFit(data        = data,
                         targetVars  = as.character(targetVars),
                         ignoreVars  = as.character(ignoreVars),
                         iterations  = as.integer(iterations),
                         sampleSizes = sampleSizes,
                         missCode    = as.integer(missCode),
                         verbose     = verbose,
                         doImp       = doImp,
                         doMcem      = doMcem,
                         doBl        = doBl)

    ## Process and check the user inputs:
    mibrrFit$processInputs()

    ## Setup the PRNG (make sure each target variable gets its own independent
    ## RNG stream):
    if(!is.null(seed))
        .lec.SetPackageSeed(rep(seed, 6))
    .lec.CreateStream(paste0("mibrrStream", c(1 : length(mibrrFit$targetVars))))
    
    ## Store Lambda's prior parameters:
    if(!doMcem) mibrrFit$setLambdaParams(l1 = as.numeric(lam1PriorPar),
                                         l2 = as.numeric(lam2PriorPar)
                                         )
    
    ## Update any user-specified control parameters:
    if(length(control) > 0) mibrrFit$setControl(control)
    
    ## Do we have any missing data:
    haveMiss <- any(mibrrFit$countMissing() > 0)

    ## Temporarily fill missing with single imputations:
    if(haveMiss) mibrrFit$simpleImpute(covsOnly = mibrrFit$fimlStarts)

    ## Are we doing any transformations?
    trans <- mibrrFit$scale || mibrrFit$center
    
    ## Compute summary statistics:
    if(trans)
        mibrrFit$computeStats(useFiml = mibrrFit$fimlStarts)
   
    if(mibrrFit$center) mibrrFit$meanCenter()
    
    ## Initialize starting values for the Gibbs sampled parameters.
    ## Important to call this before the NAs are replaced with missCode.
    mibrrFit$startParams()
    
    ## Fill remaining missing data with an integer code:
    if(mibrrFit$fimlStarts & haveMiss) mibrrFit$applyMissCode()

    mibrrFit
}# END preProcess()


## Estimate a model using MCEM:
mcem <- function(mibrrFit) {
    iters      <- mibrrFit$iterations
    totalIters <- sum(iters)
    
    for(i in 1 : totalIters) {
        if(i == 1) {
            vcat("\nBeginning MCEM 'Approximation' phase\n")
            phase <- 1
        }
        if(i == (iters[1] + 1)) {
            vcat("\nBeginning MCEM 'Tuning' phase\n")
            phase <- 2
        }
        if(i == totalIters) {
            vcat("\nSampling from the stationary posterior\n")
            phase <- 3
        }
        vcat("\n") # Beautify output
        
        ## Estimate the BEN/BL model:
        mibrrFit$doGibbs(phase)
               
        if(i < totalIters) {
            ## Print a nice message:
            if(mibrrFit$verbose) {
                check <- i > iters[1]
                cat(paste0("Doing MCEM ",
                           ifelse(check, "'Tuning'", "'Approximation'"),
                           " iteration ",
                           ifelse(check, i - iters[1], i),
                           " of ",
                           iters[as.numeric(check) + 1],
                           "\n")
                    )
            }
            
            ## Update Lambda via MCEM:
            mibrrFit$optimizeLambda(iter = i)
        }
    }# END for(i in 1 : totalIters)
    
    mibrrFit
}# END runMcem()


postProcess <- function(mibrrFit) {
    ## Uncenter the data:
    if(mibrrFit$center) mibrrFit$meanCenter(revert = TRUE)
    
    ## Replace missing values:
    if(mibrrFit$doImp) mibrrFit$applyMissCode(revert = TRUE)
    
    ## Compute the potential scale reduction factors (R-Hats) for the posterior
    ## imputation model parameters:
    if(mibrrFit$checkConv) {
        mibrrFit$computeRHats()
        mibrrFit$checkGibbsConv()
    }
    
    ## Provide some pretty names for the output objects:
    mibrrFit$nameOutput()

    ## Clean up the RNG state:
    .lec.DeleteStream(paste0("mibrrStream", c(1 : mibrrFit$nTargets)))
    
mibrrFit
}# END postProcess()
