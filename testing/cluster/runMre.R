### Title:    Run Tests to Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2019-01-01
### Modified: 2019-11-22

rm(list = ls(all = TRUE))
    
library(parallel)
library(MIBRR)
library(SURF)

resDir <- "output/"
nCores <- 2
nReps  <- 2

source("subroutines-mre.R")

## Fixed design parameters:
parms <- list()

parms$nObs    <- 100
parms$nVars   <- 10
parms$xCor    <- 0.5
parms$verbose <- TRUE
parms$doBen   <- TRUE

                                        #out <- mre(parms = parms)

## Run test in parallel:
out <- mclapply(X        = 1 : nReps,
                FUN      = mre,
                parms    = parms,
                mc.cores = nCores)

saveRDS(list(parms = parms, out = out),
        paste0(resDir,
               "mreOut_",
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S"),
               ".rds")
        )
