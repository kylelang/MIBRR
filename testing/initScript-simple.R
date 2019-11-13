### Title:    Initialize Dependencies and Parameters
### Author:   Kyle M. Lang
### Created:  2018-06-07
### Modified: 2019-11-13

library(rlecuyer)
library(MIBRR)
library(SURF)

source("simpleSubroutines.R")

## Fixed design parameters:
parms <- list()

parms$nObs       <- nObs
parms$nVars      <- nVars
parms$nTargets   <- nTargets
parms$xCor       <- xCor
parms$nReps      <- nReps
parms$streamStem <- "myStream"
parms$verbose    <- verbose
parms$mySeed     <- 235711
parms$checkKkt   <- TRUE
parms$optMeth    <- "L-BFGS-B"
parms$optBound   <- TRUE
parms$optStrict  <- TRUE
parms$iters      <- c(200, 200, 50, 50, 50, 100, 100, 1000, 1000)

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
