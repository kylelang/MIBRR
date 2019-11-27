### Title:    Run Tests to Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2019-01-01
### Modified: 2019-11-26

rm(list = ls(all = TRUE))
    
library(parallel)

verbose <- TRUE
resDir  <- "output/"
nCores  <- 2
nObs    <- 100
nVars   <- 10
pm      <- 0.0
xCor    <- 0.5
nReps   <- 8
mi      <- FALSE
sparse  <- TRUE
nPreds  <- 7

source("initScript-simple.R")

                                        #out <- testMcem(rp = 1, pm = pm, parms = parms, mi = mi)

                                        #class(out[[1]]$bl)
                                        #class(out[[2]]$bl)
                                        #class(out[[1]]$ben)
                                        #class(out[[2]]$ben)

## Run test in parallel:
time <- system.time(
    out <- mclapply(X        = 5 : nReps,
                    FUN      = testMcem,
                    pm       = pm,
                    parms    = parms,
                    nChains  = 2,
                    mi       = mi,
                    mc.cores = nCores)
)

out

lapply(out,
       function(y)
           lapply(y, function(x) {
               print(class(x$bl))
               print(class(x$ben))
           }
           )
       )

class(out[[1]]$bl)
class(out[[2]]$bl)
class(out[[1]]$ben)
class(out[[2]]$ben)

time / 60

saveRDS(list(mi = mi, parms = parms, out = out),
        paste0(resDir,
               "testMcemOut_",
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S"),
               ".rds")
        )
