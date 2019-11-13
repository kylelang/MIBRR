### Title:    Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2019-01-01
### Modified: 2019-11-13

rm(list = ls(all = TRUE))

                                        #install.packages(c("Rcpp", "RcppEigen"), repos = "http://cloud.r-project.org")

                                        #library(devtools)
                                        #install_github("kylelang/MIBRR/source/MIBRR", ref = "develop")
                                        #install_github("kylelang/SURF/source/SURF", ref = "develop")

library(parallel)

nObs        <- 100
nReps       <- 16
verbose     <- FALSE
resDir      <- "output"
clusterSize <- 2
pm          <- 0.5

source("initScript-simple.R")

## Run in serial:
                                        #miOut <- testMcem(rp = 6, pm = 0.1, parms = parms, nChains = 2)

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
                   jitterStarts = FALSE)

### Close myCluster:
stopCluster(cl)

saveRDS(list(parms = parms, out = miOut),
        paste0(resDir, "exploreMcemOut_simple_pm50_16reps_50burnInN_noJitter.rds")
        )

## Find failed reps:
check <- c()
for(i in 1 : length(miOut))
    check[i] <- class(miOut[[i]][[1]]$miben) == "try-error" |
        class(miOut[[i]][[2]]$miben) == "try-error"

which(check)
goodReps <- setdiff(1 : nReps, which(check))

### Print Stuff ###

## MCEM Chains ##

targets <- miOut[[goodReps[1]]][[1]]$miben$targetVars

## MIBEN:
par(mfcol = c(2, 2))
for(i in goodReps) {
    for(v in targets) {
        l1 <- miOut[[i]][[1]]$miben$lambdaHistory[[v]][ , 1]
        l2 <- miOut[[i]][[2]]$miben$lambdaHistory[[v]][ , 1]
        
        plot(x    = l1,
             ylim = range(l1, l2),
             type = "l",
             main = paste0("L1: ", v),
             col  = "red")
        lines(l2, col = "blue")
        
        l1 <- miOut[[i]][[1]]$miben$lambdaHistory[[v]][ , 2]
        l2 <- miOut[[i]][[2]]$miben$lambdaHistory[[v]][ , 2]
        
        plot(x    = l1,
             ylim = range(l1, l2),
             type = "l",
             main = paste0("L2: ", v),
             col  = "red")
        lines(l2, col = "blue")
    }
    readline("Hit any key to continue. ")
}

## MIBL:
par(mfcol = c(1, 2))
for(i in 1 : nReps) {
    for(v in targets) {
        l1 <- miOut[[i]][[1]]$mibl$lambdaHistory[[v]][ , 1]
        l2 <- miOut[[i]][[2]]$mibl$lambdaHistory[[v]][ , 1]
        
        plot(x    = l1,
             ylim = range(l1, l2),
             type = "l",
             main = paste0("Lambda: ", v),
             col  = "red")
        lines(l2, col = "blue")
    }
    readline("Hit any key to continue. ")
}

## Markov Chains ##

## MIBEN:
par(mfrow = c(1, 1))
for(v in with(parms, c(y, X))) {
    tmp1 <- getParams(miOut[[i]][[1]]$miben, v)
    tmp2 <- getParams(miOut[[i]][[2]]$miben, v)

    cat(paste0("\nThese are the parameters for ", v, ".\n"))
   
    plotTrace(tmp1, tmp2, "Beta")
    plotTrace(tmp1, tmp2, "Tau")
    plotTrace(tmp1, tmp2, "Sigma")
    plotTrace(tmp1, tmp2, "Lambda")
}

## MIBL:
for(v in with(parms, c(y, X))) {
    tmp1 <- getParams(miOut[[i]][[1]]$mibl, v)
    tmp2 <- getParams(miOut[[i]][[2]]$mibl, v)

    cat(paste0("\nThese are the parameters for ", v, ".\n"))
    
    plotTrace(tmp1, tmp2, "Beta")
    plotTrace(tmp1, tmp2, "Tau")
    plotTrace(tmp1, tmp2, "Sigma")
    plotTrace(tmp1, tmp2, "Lambda")
}
