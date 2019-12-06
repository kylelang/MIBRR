### Title:    Initialize Dependencies and Parameters
### Author:   Kyle M. Lang
### Created:  2018-06-07
### Modified: 2019-11-28

library(rlecuyer)
library(MIBRR)
library(SURF)

source("subroutines-simple.R")

## Fixed design parameters:
parms <- list()

parms$nObs       <- nObs
parms$nVars      <- nVars
parms$nPreds     <- nPreds
parms$xCor       <- xCor
parms$nReps      <- nReps
parms$sparse     <- sparse
parms$streamStem <- "myStream"
parms$verbose    <- verbose
parms$mySeed     <- 235711
parms$checkKkt   <- TRUE
parms$optMeth    <- "L-BFGS-B"
parms$optBound   <- TRUE
parms$optStrict  <- TRUE
parms$iters      <- c(50, 50, 10, 250, 250, 500, 500, 1000, 1000)
parms$usePcStart <- pcStart
parms$centerType <- cenType
parms$dumpPH     <- dumpPH

names(parms$iters) <- c("nMibenEmApprox",
                        "nMiblEmApprox",
                        "nEmTune",
                        "approxBurn",
                        "approxGibbs",
                        "tuneBurn",
                        "tuneGibbs",
                        "finalBurn",
                        "finalGibbs")

parms$lamStarts0 <- parms$lamStarts <- list(mibl = 7.5,  miben = c(1.0, 20.0))
