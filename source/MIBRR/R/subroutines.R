### Title:    Subroutines for the MIBRR Package
### Author:   Kyle M. Lang
### Created:  2017-11-28
### Modified: 2019-12-10

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


init <- function(penalty,
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
                 ridge,
                 verbose,
                 seed,
                 userRng,
                 nChains,
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
                         seed        = seed,
                         userRng     = userRng,
                         ridge       = ridge,
                         penalty     = as.integer(penalty),
                         nChains     = as.integer(nChains) 
                         )

    ## Process and check the user inputs:
    mibrrFit$processInputs()
    
    ## Setup the PRNG (each target variable gets an independent RNG stream):
    mibrrFit$setupRng()
    
    ## Store Lambda's prior parameters:
    if(!doMcem & penalty != 0)
        mibrrFit$setLambdaParams(l1 = as.numeric(lam1PriorPar),
                                 l2 = as.numeric(lam2PriorPar)
                                 )
    
    ## Populate the control parameters:
    if(length(control) > 0)
        MIBRR_CONTROL[names(control)] <- control

    mibrrFit$control <- MIBRR_CONTROL
        
    ## Do we have any missing data:
    haveMiss <- any(mibrrFit$missCounts > 0)

    ## Temporarily fill missing with single imputations:
    if(haveMiss) mibrrFit$simpleImpute()

    ## Initialize the 'MibrrChain' objects:
    mibrrFit$initChains()
    
    mibrrFit
}# END init()

###--------------------------------------------------------------------------###

## Estimate a model using MCEM:
mcem <- function(mibrrFit, chain) {
    iters      <- mibrrFit$iterations
    totalIters <- sum(iters)

    ## Create container for rolling parameter history:
    if(mibrrFit$control$dumpParamHistory) {
        phMem <- mibrrFit$control$phHistoryLength
        ph    <- vector("list", phMem)
    }
    
    for(i in 1 : totalIters) {
        if(i == 1) {
            vcat(paste0("\nChain ",
                        chain,
                        ": Beginning MCEM 'Approximation' phase\n")
                 )
            phase <- 1
        }
        if(i == (iters[1] + 1)) {
            vcat(paste0("\nChain ",
                        chain,
                        ": Beginning MCEM 'Tuning' phase\n")
                 )
            phase <- 2
        }
        if(i == totalIters) {
            vcat(paste0("\nChain ",
                        chain,
                        ": Sampling from the stationary posterior\n")
                 )
            phase <- 3
        }
        vcat("\n") # Beautify output
        
        ## Estimate the BEN/BL model:
        mibrrFit$chains[[chain]]$doGibbs(phase)
        
        ## Update the rolling parameter sample histories:
        if(mibrrFit$control$dumpParamHistory) 
            for (j in mibrrFit$targetVars) {
                k            <- ifelse(i %% phMem != 0, i %% phMem, phMem)
                ph[[k]][[j]] <- with(mibrrFit$gibbsOut[[j]],
                                     list(beta   = beta,
                                          tau    = tau,
                                          sigma  = sigma,
                                          lambda = mibrrFit$lambdaMat[j, ])
                                     )
                saveRDS(ph, "parameterHistory.rds")
            }
        
        if(i < totalIters) {
            ## Print a nice message:
            check <- i > iters[1]
            vcat(paste0("Chain ",
                        chain,
                        ": Doing MCEM ",
                        ifelse(check, "'Tuning'", "'Approximation'"),
                        " iteration ",
                        ifelse(check, i - iters[1], i),
                        " of ",
                        iters[as.numeric(check) + 1],
                        "\n")
                 )
                        
            ## Update Lambda via MCEM:
            optCond <- try(
                mibrrFit$chains[[chain]]$optimizeLambda(),
                silent = TRUE
            )
        }

        ## Clean up the RNG in the event of an error:
        if(class(optCond) == "try-error") {
            cleanRng(mibrrFit)
            stop(gsub("Error : ", "", optCond[[1]]), call. = FALSE)
        }
    }# END for(i in 1 : totalIters)
    
    mibrrFit
}# END mcem()

###--------------------------------------------------------------------------###

estimateModel <- function(mibrrFit) {
    for(k in 1 : mibrrFit$nChains) {
        if(mibrrFit$doMcem) mibrrFit <- mcem(mibrrFit, chain = k)
        else                mibrrFit$chains[[k]]$doGibbs()
    }
}

###--------------------------------------------------------------------------###

cleanRng <- function(mibrrFit) {
    ## Clean the RNG state:
    mibrrFit$cleanRng()
    
    ## Fix rlecuyer's random seed table when we have only 1 remaining stream:
    check <- !is.null(.lec.Random.seed.table$name) &&
        !is.matrix(.lec.Random.seed.table$Cg)
    
    if(check) .lec.Random.seed.table[1 : 4] <<-
                  lapply(.lec.Random.seed.table[1 : 4], matrix, nrow = 1)
}

###--------------------------------------------------------------------------###

postProcess <- function(mibrrFit, ...) {
    ## Extract extra arguments:
    args <- list(...)

    if(is.null(args$initOnly) || !args$initOnly) {
        ## Replace missing values:
        if(mibrrFit$doImp) mibrrFit$resetMissing()
        
        ## Compute the potential scale reduction factors (R-Hats) for the
        ## posterior imputation model parameters:
        if(mibrrFit$control$checkConv) {
            mibrrFit$computeRHats()
            mibrrFit$checkGibbsConv()
        }
        
        ## Provide some pretty names for the output objects:
        mibrrFit$nameOutput()
    }
    
    ## Clean the RNG state:
    cleanRng(mibrrFit)
        
    mibrrFit
}# END postProcess()
