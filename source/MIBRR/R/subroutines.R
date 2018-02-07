### Title:    Subroutines for the MIBRR Package
### Author:   Kyle M. Lang
### Created:  2017-NOV-28
### Modified: 2018-FEB-07

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

rm(list = ls(all = TRUE))

source("01_MibrrFit.R")
source("02_EstimationMethods.R")
source("helperFunctions.R")
source("RcppExports.R")

load("../data/mibrrExampleData.RData")

#doBl           <- FALSE
#doImp          <- TRUE
#doMcem         <- TRUE
#data           <- mibrrExampleData
#nImps          <- 100
#targetVars     <- c("y", paste0("x", c(1 : 3)))
#ignoreVars     <- "idNum"
#iterations     <- c(30, 10)
#sampleSizes    <- list(rep(25, 2),
#                      rep(250, 2),
#                      rep(500, 2)
#                      )
#missCode       <- NA
#returnConvInfo <- TRUE
#returnParams   <- FALSE
#verbose        <- TRUE
#seed           <- NULL
#control        <- list()



init <- function(doBl,
                 doImp,
                 doMcem,
                 data,
                 nImps,
                 targetVars,
                 ignoreVars,
                 iterations,
                 sampleSizes,
                 missCode,
                 verbose,
                 seed,
                 control)
{
    if(!is.null(seed)) set.seed(seed)

    ## Initialize a new MibrrFit object:
    mibrrFit <- MibrrFit(data        = data,
                         targetVars  = targetVars,
                         ignoreVars  = ignoreVars,
                         nImps       = as.integer(nImps),
                         iterations  = as.integer(iterations),
                         sampleSizes = sampleSizes,
                         missCode    = as.integer(missCode),
                         verbose     = verbose,
                         doImp       = doImp,
                         doMcem      = doMcem,
                         doBl        = doBl)

    ## Check the user inputs and resolve a set of target variables:
    mibrrFit$checkInputs()
       
    ## Update any user-specified control parameters:
    if(length(control) > 0) mibrrFit$setControl()
    
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


mibrrFit1 <- init(doBl           = FALSE,
                  doImp          = TRUE,
                  doMcem         = TRUE,
                  data           = mibrrExampleData,
                  nImps          = 100,
                  targetVars     = c("y", paste0("x", c(1 : 3))),
                  ignoreVars     = "idNum",
                  iterations     = c(30, 10),
                  sampleSizes    = list(rep(25, 2),
                                        rep(250, 2),
                                        rep(500, 2)
                                        ),
                  missCode       = NA,
                  verbose        = TRUE,
                  seed           = NULL,
                  control        = list()
                  )


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
    
mibrrFit1 <- mcem(mibrrFit1)

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
    
mibrrFit
}# END postProcess()

mibrrFit1 <- postProcess(mibrrFit1)

test <- mibrrFit1$getImpDataset()

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
