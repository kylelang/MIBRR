### Title:    Dissertation Simulation Re-Run Iteration Planning
### Author:   Kyle M. Lang
### Created:  2015-JAN-01
### Modified: 2019-MAR-05

rm(list = ls(all = TRUE))

install.packages(c("Rcpp", "RcppEigen"), repos = "http://cloud.r-project.org")

library(devtools)
install_github("kylelang/MIBRR/source/MIBRR", ref = "develop")
                                        #install_github("kylelang/SURF/source/SURF", ref = "develop")

library(parallel)

expNum      <- 1

sparse      <- FALSE
nImps       <- 100
nObs        <- 100
startRep    <- 1
stopRep     <- 12
clusterSize <- 3
outDir      <- NULL

resDir <- "/home/kmlang/sg/research/active/miben/dissRerun/results/isbaRerun/iterPlan/"

source("initScript_working.R")


## Create the cluster:
cl <- makeCluster(clusterSize)

clusterCall(cl = cl, fun = source, file = "subroutines.R")
clusterCall(cl      = cl,
            fun     = applyLib,
            pkgList = c("rlecuyer", "mice", "mitools", "MIBRR", "SURF")
            )

### Apply doMiOnly() in parallel:
miOut <- parLapply(cl      = cl,
                   X       = startRep : stopRep,
                   fun     = iterPlan,
                   pm      = 0.1,
                   control = control,
                   parms   = parms,
                   doMice  = FALSE)

### Close myCluster:
stopCluster(cl)

saveRDS(list(control = control, out = miOut),
        paste0(resDir,
               "iterPlanningOut_", ifelse(control$sparse, "sparse", "dense"),
               "_exp",             expNum,
               ".rds")
        )

tmp <- readRDS(paste0(resDir,
                      "iterPlanningOut_", ifelse(control$sparse, "sparse", "dense"),
                      "_exp",             expNum,
                      ".rds")
               )

miOut <- tmp$out

## Find failed reps:
check <- c()
for(i in 1 : length(miOut))
    check[i] <- class(miOut[[i]][[1]]$miben) == "try-error" |
        class(miOut[[i]][[2]]$miben) == "try-error"

which(check)

### Print Stuff ###

## Markov Chains ##

## MIBEN:
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


## MCEM Chains ##

i <- 9
v <- "x2"

par(mfrow = c(1, 1))

l1 <- miOut[[i]][[1]]$miben$lambdaHistory[[v]][ , 2]
l2 <- miOut[[i]][[2]]$miben$lambdaHistory[[v]][ , 2]

plot(x    = l1,
     ylim = range(l1, l2),
     type = "l",
     ylab = "Lambda2",
     xlab = "Iteration",
     col  = "red")
lines(l2, col = "blue")


## MIBEN:
par(mfcol = c(2, 3))
for(v in with(parms, c(y, X))) {
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

## MIBL:
par(mfcol = c(1, 3))
for(v in with(parms, c(y, X))) {
    l1 <- miOut[[i]][[1]]$mibl$lambdaHistory[[v]][ , 1]
    l2 <- miOut[[i]][[2]]$mibl$lambdaHistory[[v]][ , 1]
    
    plot(x    = l1,
         ylim = range(l1, l2),
         type = "l",
         main = paste0("Lambda: ", v),
         col  = "red")
    lines(l2, col = "blue")
}

plot(miOut[[i]]$naiveMice, layout = c(2, 3))
plot(miOut[[i]]$quickMice, layout = c(2, 3))
plot(miOut[[i]]$trueMice,  layout = c(2, 3))

saveRDS(list(dense = miOutD, sparse = NA),
        file = paste0(outDir, "exp3_iter_planning_results.rds")
        )
miOut <- readRDS(paste0(outDir, "exp2_iter_planning_results.rds"))

miblS <- miblD <- mibenS <- mibenD <- list()
for(i in 1 : 4) {
    miblS[[i]]  <- getParams(miOut$sparse[[i]]$mibl,  "y")$lambda
    miblD[[i]]  <- getParams(miOut$dense[[i]]$mibl,   "y")$lambda
    mibenS[[i]] <- getParams(miOut$sparse[[i]]$miben, "y")$lambda
    mibenD[[i]] <- getParams(miOut$dense[[i]]$miben,  "y")$lambda
}

miblS
miblD
mibenS
mibenD

## Time a single replication:
repTime <- system.time(
    mclapply(X        = c(startRep : stopRep),
             FUN      = goBabyGo,
             control  = control,
             parms    = parms,
             mc.cores = clusterSize)
)

saveRDS(repTime, paste0(outDir, "exp3_dense_rep_time.rds"))
                                        #repTime <- readRDS(paste0(outDir, "exp2_sparse_rep_time.rds"))

getRunTime(repTime[["elapsed"]], 500, 60)

## Check results:
nA <- 120
pm <- 30
sp <- "dense"
rp <- 2

resList <- readRDS(
    paste0(outDir, "simRes_", sp, "_v", nA, "_pm", pm, "_rep", rp, ".rds")
)

resList$comp
resList$naiveMice
resList$quickMice
resList$trueMice
resList$miben
resList$mibl
