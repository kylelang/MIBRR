### Title:    Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2019-01-01
### Modified: 2019-11-20

rm(list = ls(all = TRUE))

                                        #install.packages(c("Rcpp", "RcppEigen"), repos = "http://cloud.r-project.org")

                                        #library(devtools)
                                        #install_github("kylelang/MIBRR/source/MIBRR", ref = "develop")
                                        #install_github("kylelang/SURF/source/SURF", ref = "develop")

library(parallel)

nReps       <- 16
verbose     <- FALSE
resDir      <- "output"
clusterSize <- 2
nObs        <- 100
nVars       <- 10
nTargets    <- 4
pm          <- 0.0
xCor        <- 0.5

source("initScript-simple.R")

## Run in serial:
                                        #miOut <- testMcem(rp = 1, pm = 0.1, parms = parms, nChains = 2, mi = FALSE)

## Create the cluster:
cl <- makeCluster(clusterSize)

clusterCall(cl = cl, fun = source, file = "simpleSubroutines.R")
clusterCall(cl      = cl,
            fun     = applyLib,
            pkgList = c("rlecuyer", "MIBRR", "SURF")
            )

### Apply iterPlan() in parallel:
miOut <- parLapply(cl           = cl,
                   X            = 1 : nReps,
                   fun          = testMcem,
                   pm           = pm,
                   parms        = parms,
                   jitterStarts = TRUE,
                   mi           = FALSE)

### Close myCluster:
stopCluster(cl)

saveRDS(list(parms = parms, out = miOut),
        paste0(resDir, "exploreMcemOut_simple_ben_p10_xCor0_500gibbsN.rds")
        )

                                        #miOut <- readRDS(paste0(resDir, "exploreMcemOut_simple_ben_p10_xCor0.rds"))$out

## Find failed reps:
check <- c()
for(i in 1 : length(miOut))
    check[i] <- class(miOut[[i]][[1]]$ben) == "try-error" |
        class(miOut[[i]][[2]]$ben) == "try-error"

which(check)
goodReps <- setdiff(1 : nReps, which(check))

### Print Stuff ###

## MCEM Chains ##

## MIBEN:
par(mfcol = c(1, 2))
for(i in goodReps) {
    l1 <- miOut[[i]][[1]]$ben$lambdaHistory[[1]][ , 1]
    l2 <- miOut[[i]][[2]]$ben$lambdaHistory[[1]][ , 1]
    
    plot(x    = l1,
         ylim = range(l1, l2),
         type = "l",
         main = "L1",
         col  = "red")
    lines(l2, col = "blue")
    
    l1 <- miOut[[i]][[1]]$ben$lambdaHistory[[1]][ , 2]
    l2 <- miOut[[i]][[2]]$ben$lambdaHistory[[1]][ , 2]
    
    plot(x    = l1,
         ylim = range(l1, l2),
         type = "l",
         main = "L2",
         col  = "red")
    lines(l2, col = "blue")
    
    readline("Hit any key to continue. ")
}

## MIBL:
par(mfcol = c(1, 1))
for(i in 1 : nReps) {
    l1 <- miOut[[i]][[1]]$bl$lambdaHistory[[1]][ , 1]
    l2 <- miOut[[i]][[2]]$bl$lambdaHistory[[1]][ , 1]
    
    plot(x    = l1,
         ylim = range(l1, l2),
         type = "l",
         main = "Lambda",
         col  = "red")
    lines(l2, col = "blue")

    readline("Hit any key to continue. ")
}

## Markov Chains ##

length(miOut)

## MIBEN:
par(mfrow = c(1, 1))
for(i in goodReps) {
    tmp1 <- getParams(miOut[[i]][[1]]$ben, 1)
    tmp2 <- getParams(miOut[[i]][[2]]$ben, 1)
    
    plotTrace(tmp1, tmp2, "Beta")
    plotTrace(tmp1, tmp2, "Tau")
    plotTrace(tmp1, tmp2, "Sigma")
    
    readline("Hit any key to continue. ")
}

## MIBL:
for(i in 1 : nReps) {
    tmp1 <- getParams(miOut[[i]][[1]]$bl, 1)
    tmp2 <- getParams(miOut[[i]][[2]]$bl, 1)
    
    plotTrace(tmp1, tmp2, "Beta")
    plotTrace(tmp1, tmp2, "Tau")
    plotTrace(tmp1, tmp2, "Sigma")

    readline("Hit any key to continue. ")
}
