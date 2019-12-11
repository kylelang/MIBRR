### Title:    Test Reference Classes for MIBRR Package
### Author:   Kyle M. Lang
### Created:  2019-10-09
### Modified: 2019-12-11

rm(list = ls(all = TRUE))

source("../source/MIBRR/R/00_MibrrSamples.R")
source("../source/MIBRR/R/01_MibrrChain.R")
source("../source/MIBRR/R/02_MibrrFit.R")

source("../source/MIBRR/R/initControl.R")
source("../source/MIBRR/R/helperFunctions.R")
source("../source/MIBRR/R/subroutines.R")

library(mvtnorm)
library(MIBRR)
library(rlecuyer)
library(optimx)

source("../source/MIBRR/R/exportedPrimaryFunctions.R")

                                        #library(devtools)
                                        #install_github("kylelang/MIBRR/source/MIBRR", ref = "develop")

data(mibrrExampleData)

out <- miben(data        = mibrrExampleData,
             targetVars  = c("y", paste0("x", c(1 : 3))),
             ignoreVars  = "idNum",
             iterations  = c(30L, 10L),
             sampleSizes = list(rep(50, 2), rep(100, 2), rep(500, 2)),
             doMcem      = TRUE,
             verbose     = TRUE,
             seed        = 235711,
             nChains     = 2L,
             nCores      = 1L,
             userRng     = "",
             control     = list(checkConv = TRUE)
             )

length(out$chains)

out$chains[[1]]
out$chains[[2]]

out$getImpDataset()

mibrrFit <- init(data         = mibrrExampleData,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                 iterations   = c(30L, 10L),
                 sampleSizes  = list(rep(50, 2), rep(100, 2), rep(500, 2)),
                 doImp        = TRUE,
                 doMcem       = TRUE,
                 verbose      = TRUE,
                 seed         = 235711,
                 penalty      = 2L,
                 nChains      = 3L,
                 missCode     = NA,
                 lam1PriorPar = NA,
                 lam2PriorPar = NA,
                 ridge        = 0.0,
                 userRng      = "",
                 control      = list(checkConv = TRUE)
                 )

for(k in 1 : 3)
    mibrrFit <- mcem(mibrrFit, chain = k)
mibrrFit <- postProcess(mibrrFit)

mibrrFit$rHats

imps <- list()
for(m in 1 : 5)
    imps[[m]] <- mibrrFit$getImpDataset()

all.equal(imps[[1]], imps[[2]])
all.equal(imps[[1]], imps[[3]])
all.equal(imps[[1]], imps[[4]])
all.equal(imps[[1]], imps[[5]])

all.equal(imps[[2]], imps[[3]])
all.equal(imps[[2]], imps[[4]])
all.equal(imps[[2]], imps[[5]])

all.equal(imps[[3]], imps[[4]])
all.equal(imps[[3]], imps[[5]])

all.equal(imps[[4]], imps[[5]])


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


tmp <- list(matrix(rnorm(100), 20, 5), matrix(rnorm(100), 20, 5))

do.call(rbind, tmp)
