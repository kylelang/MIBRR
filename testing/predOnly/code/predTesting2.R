### Title:    Test MIBEN/MIBL Prediction
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2016-MAY-16

rm(list = ls(all = TRUE))

outDir <- "../output/"

library(mice)
library(mitools)
library(mvtnorm)
library(mibrr)
library(parallel)
library(Matrix)

## Test BEN and BL:
alpha  <- 50
nPreds <- 90
nReps  <- 500

condMat <- expand.grid(center = c(TRUE, FALSE), scale = c(TRUE, FALSE))

means  <- runif(nPreds, 0, 100)
scales <- runif(nPreds, 1, 100)
collin <- 0.3

outMeans <- matrix(NA, 4, 3)
outSds <- matrix(NA, 4, 3)
outMse <- list()
for(i in 1 : nrow(condMat)) {
    print(paste0("Doing condition set ", i))
    print(" ")
    
    parms <- list()
    parms$nObs        <- 100
    parms$nPreds      <- nPreds
    parms$r2          <- 0.5
    parms$collin      <- collin
    parms$beta        <- matrix(
        c(alpha, runif(nPreds / 2, 0.3, 0.6), rep(0, nPreds / 2))
    )
                                        #parms$beta        <- matrix(c(alpha, runif(nPreds, 0.3, 0.6)))
    parms$means       <- means
    parms$scales      <- scales
    parms$center      <- condMat[i, "center"]
    parms$scale       <- condMat[i, "scale"]
    parms$useClassic  <- FALSE
    parms$simpleInt   <- FALSE
    parms$verbose     <- FALSE
    parms$postN       <- 1000
    parms$approxIters <- 200
    parms$approxN     <- 100
    parms$tuneIters   <- 20
    parms$tuneN       <- 250

    mseList <- mclapply(c(1 : nReps),
                        FUN = testFun,
                        parms = parms,
                        mc.cores = 4)
  
    mseMat <- do.call(rbind, mseList)
    outMse[[i]] <- mseMat
    outMeans[i, ] <- colMeans(mseMat)
    outSds[i, ] <- apply(mseMat, 2, sd)
}


outList <- list(
    mean = cbind(condMat, outMeans),
    sd =  cbind(condMat, outSds)
)

                                        #outList$mean[ , c(3, 4)] / outList$mean[ , 5]
                                        #outList$sd[ , c(3, 4)] / outList$sd[ , 5]

saveRDS(outList,
        file = paste0(outDir, "sparse_90pred_outList1.rds")
        )
