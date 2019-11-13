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

                                        #parms$nImps      <- nImps
parms$nObs       <- nObs
parms$cx <- 0.3
                                        #parms$y          <- "y"
                                        #parms$X          <- paste0("x", 1 : 2)
                                        #parms$beta       <- rep(1.0, length(parms$X) + 1)
                                        #parms$zeta       <- 0.3
                                        #parms$r2         <- list(X = 0.25, y = 0.5)
                                        #parms$muX        <- rep(1.0, length(parms$X))
                                        #parms$pmVec      <- c(0.1, 0.2, 0.3)
                                        #parms$marSnr     <- 5.0
                                        #parms$marPattern <- c(y = "low", x1 = "high", x2 = "low")
                                        #parms$minPredCor <- 0.3
                                        #parms$pmAux      <- 0.0
                                        #parms$outDir     <- outDir
parms$nReps      <- nReps
parms$streamStem <- "myStream"
parms$verbose    <- verbose
parms$mySeed     <- 235711
                                        #parms$saveParams <- TRUE
parms$checkKkt   <- TRUE
parms$optMeth    <- "L-BFGS-B"
parms$optBound   <- TRUE
                                        #parms$miceMeth   <- "norm"
parms$optStrict  <- TRUE
                                        #parms$mcem       <- TRUE

## Variable design parameters:
                                        #control <- list()

                                        #control$expNum    <- expNum
                                        #control$sparse    <- sparse
                                        #control$nAux      <- switch(expNum, 10, 60, 120)
                                        #control$marPreds  <- paste0("z", 1 : (control$nAux / 5))
                                        #control$miceIters <- switch(expNum, 10, 10, 10)
                                        #control$miceRidge <- switch(expNum, 1e-5, 1e-1, 1e-1)
                                        #control$vanRidge  <- 1e-2
                                        #control$l1Pars    <- c(1.5, 0.05)
                                        #control$l2Pars    <- c(15.0, 0.025)
parms$iters     <- c(500, 500, 100, 25, 25, 100, 100, 1000, 1000)

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
