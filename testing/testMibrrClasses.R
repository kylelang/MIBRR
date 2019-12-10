### Title:    Test Reference Classes for MIBRR Package
### Author:   Kyle M. Lang
### Created:  2019-10-09
### Modified: 2019-12-10

rm(list = ls(all = TRUE))

source("../source/MIBRR/R/00_MibrrSamples.R")
source("../source/MIBRR/R/01_MibrrChain.R")
source("../source/MIBRR/R/02_MibrrFit.R")

source("../source/MIBRR/R/initControl.R")

library(mvtnorm)
library(MIBRR)
library(rlecuyer)
library(optimx)
                                        #library(devtools)
                                        #install_github("kylelang/MIBRR/source/MIBRR", ref = "develop")

data(mibrrExampleData)

## Estimate the mode of a continous vector:
numMode <- function(x) {
    dens <- density(x, na.rm = TRUE)
    dens$x[which.max(dens$y)]
}

cast <- function(obj, type)
    eval(call(paste0("as.", type), obj))

setControl <- function(x, where = .self) {
    ## Get the fields for the current class:
    fields <- getRefClass(class(where))$fields()
    
    ## Assign the control list entries to the correct classes:
    for(n in names(x))
        if(n %in% names(fields))
            where$field(n, cast(x[n], fields[n]))
}

mibrrFit <- MibrrFit(data        = mibrrExampleData,
                     targetVars  = c("y", paste0("x", c(1 : 3))),
                     ignoreVars  = "idNum",
                     iterations  = c(30L, 10L),
                     sampleSizes = list(rep(50, 2), rep(100, 2), rep(500, 2)),
                     doImp       = TRUE,
                     doMcem      = TRUE,
                     verbose     = TRUE,
                     seed        = 235711,
                     penalty     = 2L,
                     nChains     = 2L)

## Process and check the user inputs:
mibrrFit$processInputs()

## Setup the PRNG (each target variable gets an independent RNG stream):
mibrrFit$setupRng()

mibrrFit$control <- MIBRR_CONTROL

## Do we have any missing data:
haveMiss <- any(mibrrFit$missCounts > 0)

## Temporarily fill missing with single imputations:
if(haveMiss) mibrrFit$simpleImpute()

## Initialize the 'MibrrChain' objects:
mibrrFit$initChains()

mibrrFit$chains[[1]]$doGibbs()
mibrrFit$chains[[1]]$optimizeLambda()
mibrrFit$chains[[1]]$parameters[["y"]]

fit$fields
fit$methods()

source("../source/MIBRR/R/control0.R")

fitFields     <- getRefClass("MibrrFit")$fields()
chainFields   <- getRefClass("MibrrChain")$fields()
samplesFields <- getRefClass("MibrrSamples")$fields()

?expression
?call

tmp <- call("as.matrix", 2)
eval(tmp)

type <- "matrix"
eval(call(paste0("as.", type), rep(3, 3)))

eval(tmp)

tmp(rep(3, 10))


    where$field(
    
names(control) %in% names(tmp)

ls(fit)

fit$chains[[1]]$doGibbs()
chain$optimizeLambda()

chain$parameters[["y"]]$lambda1

## MCEM estimation:
                                        #mibenOut <- miben(data       = mibrrExampleData,
                                        #                  iterations = c(30, 10),
                                        #                  targetVars = c("y", paste0("x", c(1 : 3))),
                                        #                  ignoreVars = "idNum")

dat1 <- mibrrExampleData[setdiff(colnames(mibrrExampleData), "idNum")]



missList <- lapply(dat1, function(x) which(is.na(x)) - 1)

dat1[is.na(dat1)] <- 0.0
                                        #missList

                                        #nrow(dat1) - sapply(missList, length)

nChains <- 2
targetVars <- c("y", paste0("x", c(1 : 3)))
nTargets <- length(targetVars)

## Generate the l'ecuyer RNG streams:
streams <<- c("mibrrStream0",
              paste0("c",
                     rep(1 : nChains, each = nTargets),
                     targetVars)
              )

.lec.CreateStream(streams)

## Set the 'master' RNG stream:
rng0 <<- .lec.CurrentStream("mibrrStream0")

chain <- MibrrChain(chain       = 1L,
                    data        = dat1,
                    targetVars  = c("y", paste0("x", c(1 : 3))),
                    iterations  = c(30L, 10L),
                    sampleSizes = list(rep(50, 2), rep(100, 2), rep(500, 2)),
                    missList    = missList,
                    doMcem      = TRUE,
                    verbose     = TRUE,
                    centerType  = "mode",
                    penalty     = 2L)

chain$doGibbs()
chain$optimizeLambda()

chain$parameters[["y"]]$lambda1

    sams <- MibrrSamples(target      = "y",
                     predVars    = paste0("x", 1 : 10),
                     iterations  = 50L,
                     targetScale = 0.25,
                     doMcem      = TRUE,
                     penalty     = 2L)

sams

sams$incIter()

sams$iter

sams$setLambdas(c(2, 5))

sams$lambda1
sams$lambda2

sams$target

sams$setSamples(mibenOut$gibbsOut)

sams$sigma
sams$imps

names(mibenOut$gibbsOut[["y"]])

length(mibenOut$gibbsOut[["y"]][["ppSams"]])
