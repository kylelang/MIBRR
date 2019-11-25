### Title:    Initialize Dependencies and Parameters
### Author:   Kyle M. Lang
### Created:  2018-06-07
### Modified: 2019-11-25

library(MIBRR)
library(SURF)

source("subroutines-mre.R")

## Fixed design parameters:
parms <- list()

parms$nObs    <- nObs
parms$nVars   <- nVars
parms$xCor    <- xCor
parms$verbose <- verbose
parms$doBen   <- TRUE
