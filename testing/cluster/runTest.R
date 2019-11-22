### Title:    Run Tests to Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2019-01-01
### Modified: 2019-11-22

rm(list = ls(all = TRUE))
    
library(parallel)

verbose <- FALSE
resDir  <- "output/"
nCores  <- 2
nObs    <- 100
nVars   <- 10
pm      <- 0.0
xCor    <- 0.5
nReps   <- 2
mi      <- FALSE

source("initScript-simple.R")

## Run test in parallel:
time <- system.time(
    out <- mclapply(X        = 1 : nReps,
                    FUN      = testMcem,
                    pm       = 0.0,
                    parms    = parms,
                    nChains  = 2,
                    mi       = mi,
                    mc.cores = nCores)
)

saveRDS(list(mi = mi, parms = parms, out = out),
        paste0(resDir,
               "testMcemOut_",
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S"),
               ".rds")
        )

time / 60
