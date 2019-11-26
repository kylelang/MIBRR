### Title:    Run Tests to Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2019-01-01
### Modified: 2019-11-26

rm(list = ls(all = TRUE))
    
library(parallel)

verbose <- FALSE
resDir  <- "output/"
nCores  <- 2
nObs    <- 1000
nVars   <- 10
pm      <- 0.0
xCor    <- 0.5
nReps   <- 2
mi      <- FALSE
sparse  <- TRUE
nPreds  <- 5

source("initScript-simple.R")

out <- testMcem(rp = 1, pm = pm, parms = parms, mi = mi)

## Run test in parallel:
time <- system.time(
    out <- mclapply(X        = 1 : nReps,
                    FUN      = testMcem,
                    pm       = pm,
                    parms    = parms,
                    nChains  = 2,
                    mi       = mi,
                    mc.cores = nCores)
)

time / 60

saveRDS(list(mi = mi, parms = parms, out = out),
        paste0(resDir,
               "testMcemOut_",
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S"),
               ".rds")
        )
