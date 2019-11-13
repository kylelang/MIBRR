### Title:    Simple Testing Subroutines
### Author:   Kyle M. Lang
### Created:  2019-11-13
### Modified: 2019-11-13
### Purpose:  This file contains the subroutines to help debug the MCEM
###           estimation in MIBRR

###--------------------------------------------------------------------------###

                                        #library(SURF)
                                        #?simCovData
                                        #?imposeMissData

## Simulate some very simple data:
genSimpleData <- function(n, cov, pm) {
    ## Generate complete data:
    dat0 <- SURF::simCovData(nObs = n, sigma = cov, nVars = 3)

    ## Impose missing:
    r1 <- dat0$x3 < quantile(dat0$x3, pm)
    r2 <- dat0$x3 > quantile(dat0$x3, (1 - pm))
    
    dat0[r1, "x1"] <- NA
    dat0[r2, "x2"] <- NA

    dat0
}

###--------------------------------------------------------------------------###

## Estimate the imputation models:
runMibrr <- function(data, parms) {

    ## Define imputation targets:
    targets <- colnames(data)[colMeans(is.na(data)) > 0]
    
    ## Advance the RNG stream:
    .lec.ResetNextSubstream(parms$rngStream)
    
    ## Run MIBEN:
    mibenOut <- try(
        miben(data        = data,
              targetVars  = targets,
              iterations  = as.numeric(
                  parms$iters[c("nMibenEmApprox", "nEmTune")]
              ),
              sampleSizes = list(
                  as.numeric(parms$iters[c("approxBurn", "approxGibbs")]),
                  as.numeric(parms$iters[c("tuneBurn", "tuneGibbs")]),
                  as.numeric(parms$iters[c("finalBurn", "finalGibbs")])
              ),
              seed        = parms$rngStream,
              userRng     = parms$rngStream,
              verbose     = parms$verbose,
              control     = list(optCheckKkt    = parms$checkKkt,
                                 optMethod      = parms$optMeth,
                                 optBoundLambda = parms$optBound,
                                 optStrict      = parms$optStrict,
                                 lambda1Starts  = parms$lamStarts$miben[[1]],
                                 lambda2Starts  = parms$lamStarts$miben[[2]])
              )
    )
    
    ## Advance the RNG stream:
    .lec.ResetNextSubstream(parms$rngStream)

    ## Run MIBL:
    miblOut <- try(
        mibl(data        = data,
             targetVars  = targets,
             iterations  = as.numeric(
                 parms$iters[c("nMiblEmApprox", "nEmTune")]
             ),
             sampleSizes = list(
                 as.numeric(parms$iters[c("approxBurn", "approxGibbs")]),
                 as.numeric(parms$iters[c("tuneBurn", "tuneGibbs")]),
                 as.numeric(parms$iters[c("finalBurn", "finalGibbs")])
             ),
             seed        = parms$rngStream,
             userRng     = parms$rngStream,
             verbose     = parms$verbose,
             control     = list(lambda1Starts = parms$lamStarts$mibl)
             )
    )
    
    ## Return a combined output object:
    list(miben = mibenOut, mibl = miblOut)
}

###--------------------------------------------------------------------------###

## Only do MI for a simple level of PM. Used for iteration planning.
testMcem <- function(rp, pm, parms, nChains = 2) {
    ## Store the current rep index:
    parms$rp <- rp
    
    ## Initialize the L'Ecuyer random number streams:
    .lec.SetPackageSeed(rep(parms$mySeed, 6))
    .lec.CreateStream(paste0(parms$streamStem, c(1 : parms$nReps)))
    
    ## Associate this replication with it's own stream:
    parms$rngStream <- paste0(parms$streamStem, rp)
    .lec.CurrentStream(parms$rngStream)
    
    ## Simulate the complete data:
    dat0 <- genSimpleData(n = parms$nObs, cov = parms$cx, pm = pm)
    
    out <- list()
    for(i in 1 : nChains) {
        ## Perturb starting values:
        parms$lamStarts <- lapply(parms$lamStarts0,
                                  function(x) x + rnorm(length(x), 0, 0.5 * x)
                                  )
        
        ## Apply over percents missing:
        out[[i]] <- runMibrr(data = dat0, parms = parms)
    }
    
    ## Clean up the RNG state:
    .lec.CurrentStreamEnd()
    .lec.DeleteStream(paste0(parms$streamStem, 1 : parms$nReps))
    
    out
}

###--------------------------------------------------------------------------###
