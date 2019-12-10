source("00_MibrrSamples.R")
source("01_MibrrChain.R")

library(mvtnorm)
library(MIBRR)
library(devtools)

install_github("kylelang/MIBRR/source/MIBRR", ref = "develop")

?miben

data(mibrrExampleData)

## MCEM estimation:
mibenOut <- miben(data       = mibrrExampleData,
                  iterations = c(30, 10),
                  targetVars = c("y", paste0("x", c(1 : 3))),
                  ignoreVars = "idNum")

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
