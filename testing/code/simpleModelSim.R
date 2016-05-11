## Title:    Simulation using a simple regression model
## Author:   Kyle M. Lang
## Created:  2016-MAY-01
## Modified: 2016-MAY-04
## Purpose:  Run a miben simulation using a very simple regression model with
##           absolutely no need for regularization, to address reviewer comments

rm(list = ls(all = TRUE))
                                        #softDir <- "../software/"
                                        #install.packages(paste0(softDir, "mibrr_0.0.0.tar.gz"),
                                        #                 repos = NULL,
                                        #                 type = "source")

library(mibrr)
library(mvtnorm)
library(mitools)
library(mice)
library(parallel)
source("simpleSimFunctions.R")



X <- as.matrix(dat1[ , c("sysRac", "indRac", "merit")])

aMat <- crossprod(X)

test1 <- aMat * 3
test2 <- aMat %*% diag(rep(3, 3))

test1 - test2


nReps <- 200
cSize <- 8

parms <- list()
parms$nObs       <- 500
parms$pm         <- 0.30
parms$auxVars    <- c("indRac", "merit")
parms$incompVars <- c("policy", "sysRac", "polAffil", "revDisc")
parms$marType    <- c("lower", "upper", "tails", "center")
parms$testForm   <- as.formula("policy ~ polAffil + sysRac")
parms$dataDir    <- "../data/"
parms$resDir     <- "../results4/"

dat1 <- cleanData2(parms)
parms$sigma <- cov(dat1)
parms$mu <- colMeans(dat1)

### Run the simulation:
outList <- mclapply(c(1 : nReps), FUN = doRep, parms = parms, mc.cores = cSize)

## Summarize the results:
coefList <- seList <- list(comp  = matrix(NA, nReps, 3),
                           miben = matrix(NA, nReps, 3),
                           mice  = matrix(NA, nReps, 3)
                           )
failCount <- 0
for(rp in 1 : nReps) {
    tmp <- readRDS(paste0(parms$resDir, "outFile", rp, ".rds"))
    if(!is.null(tmp) & class(tmp) != "character") {
        coefList$comp[rp, ]  <- coef(tmp$comp)
        coefList$miben[rp, ] <- coef(tmp$miben)
        coefList$mice[rp, ]  <- coef(tmp$mice)
        
        seList$comp[rp, ]  <- sqrt(diag(vcov(tmp$comp)))
        seList$miben[rp, ] <- sqrt(diag(vcov(tmp$miben)))
        seList$mice[rp, ]  <- sqrt(diag(vcov(tmp$mice)))
    } else {
        failCount <- failCount + 1
    }
}

failCount

mibenCoefs <- coefList$miben
miceCoefs <- coefList$mice
compCoefs <- coefList$comp

mibenMeans <- colMeans(mibenCoefs, na.rm = TRUE)
miceMeans <- colMeans(miceCoefs, na.rm = TRUE)
compMeans <- colMeans(compCoefs, na.rm = TRUE)

prbMiben <- 100 * (mibenMeans - compMeans) / abs(compMeans)
prbMice <- 100 * (miceMeans - compMeans) / abs(compMeans)

prbMiben
prbMice

sbMiben <-
    (colMeans(mibenCoefs) - colMeans(compCoefs)) /
        apply(mibenCoefs, 2, sd)
sbMice <-
    (colMeans(miceCoefs) - colMeans(compCoefs)) /
        apply(miceCoefs, 2, sd)

sbMiben
sbMice
