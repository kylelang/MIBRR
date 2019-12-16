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
source("../source/MIBRR/R/exportedHelperFunctions.R")

library(mvtnorm)
library(MIBRR)
library(rlecuyer)
library(optimx)

source("../source/MIBRR/R/exportedPrimaryFunctions.R")

                                        #library(devtools)
                                        #install_github("kylelang/MIBRR/source/MIBRR", ref = "develop")

data(mibrrExampleData)

mibenOut <- miben(data        = mibrrExampleData,
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

mibenOut$chains[[1]]$parameters[["y"]]$logLik

pars <- getParams(mibenOut, "y", mix = TRUE)
pars$beta

ls(mibenOut)

lams <- getLambda(mibenOut, "y")

plot(lams[[1]]$conv$ll, type = "l")
lines(lams[[2]]$conv$ll, col = "red")

MIBRR:::plotLambda(mibenOut, "x2", TRUE)

?miben
?ben

mibenOut <- miben(data         = mibrrExampleData,
                  targetVars   = c("y", paste0("x", c(1 : 3))),
                  ignoreVars   = "idNum",
                  sampleSizes  = c(500, 500),
                  doMcem       = FALSE,
                  lam1PriorPar = c(1.0, 0.1),
                  lam2PriorPar = c(1.0, 0.1)
                  )
mibenOut$rHats

out <- miben(data         = mibrrExampleData,
             targetVars   = c("y", paste0("x", c(1 : 3))),
             ignoreVars   = "idNum",
             sampleSizes  = rep(500, 2),
             lam1PriorPar = c(1.0, 0.1),
             lam2PriorPar = c(1.0, 0.1),
             doMcem       = FALSE,
             verbose      = TRUE,
             seed         = 235711,
             nChains      = 2L,
             nCores       = 1L,
             userRng      = "",
             control      = list(checkConv = TRUE)
             )

out <- ben(data        = mibrrExampleData,
           y           = "y",
           X           = paste0("x", c(1 : 3)),
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

out <- mibl(data        = mibrrExampleData,
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

miblOut <- mibl(data         = mibrrExampleData,
                targetVars   = c("y", paste0("x", c(1 : 3))),
                ignoreVars   = "idNum",
                sampleSizes  = c(500, 500),
                doMcem       = FALSE,
                lam1PriorPar = c(1.0, 0.1)
                )
pars <- getParams(miblOut, "y")
pars

out <- mibl(data        = mibrrExampleData,
            targetVars  = c("y", paste0("x", c(1 : 3))),
            ignoreVars  = "idNum",
            iterations  = c(30L, 10L),
            sampleSizes = list(rep(50, 2), rep(100, 2), rep(500, 2)),
            doMcem      = FALSE,
            verbose     = TRUE,
            seed        = 235711,
            nChains     = 2L,
            nCores      = 1L,
            userRng     = "",
            control     = list(checkConv = TRUE)
            )

out <- bl(data        = mibrrExampleData,
          y           = "y",
          X           = paste0("x", c(1 : 3)),
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

out <- vanilla(data        = mibrrExampleData,
               targetVars  = c("y", paste0("x", c(1 : 3))),
               ignoreVars  = "idNum",
               sampleSizes = rep(500, 2),
               verbose     = TRUE,
               seed        = 235711,
               nChains     = 2L,
               nCores      = 1L,
               userRng     = "",
               control     = list(checkConv = TRUE)
               )

pars <- getParams(out, "y")
pars

out$chains[[1]]$parameters[["y"]]$beta

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

?optimx

## Show multiple outputs of optimx using all.methods
                                        # genrose function code
genrose.f<- function(x, gs=NULL){ # objective function
    ## One generalization of the Rosenbrock banana valley function (n parameters)
    n <- length(x)
    if(is.null(gs)) { gs=100.0 }
    fval<-1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2)
    return(fval)
}

genrose.g <- function(x, gs=NULL){
                                        # vectorized gradient for genrose.f
                                        # Ravi Varadhan 2009-04-03
    n <- length(x)
    if(is.null(gs)) { gs=100.0 }
    gg <- as.vector(rep(0, n))
    tn <- 2:n
    tn1 <- tn - 1
    z1 <- x[tn] - x[tn1]^2
    z2 <- 1 - x[tn]
    gg[tn] <- 2 * (gs * z1 - z2)
    gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
    return(gg)
}

genrose.h <- function(x, gs=NULL) { ## compute Hessian
    if(is.null(gs)) { gs=100.0 }
    n <- length(x)
    hh<-matrix(rep(0, n*n),n,n)
    for (i in 2:n) {
        z1<-x[i]-x[i-1]*x[i-1]
        z2<-1.0-x[i]
        hh[i,i]<-hh[i,i]+2.0*(gs+1.0)
        hh[i-1,i-1]<-hh[i-1,i-1]-4.0*gs*z1-4.0*gs*x[i-1]*(-2.0*x[i-1])
        hh[i,i-1]<-hh[i,i-1]-4.0*gs*x[i-1]
        hh[i-1,i]<-hh[i-1,i]-4.0*gs*x[i-1]
    }
    return(hh)
}

library(optimx)

startx<-4*seq(1:10)/3.
ans8<-optimx(startx,fn=genrose.f,gr=genrose.g, hess=genrose.h, method = "BFGS", gs=10)

class(ans8)

ans8$value

ans8[, "gevals"]
ans8["spg", ]
summary(ans8, par.select = 1:3)
summary(ans8, order = value)[1, ] # show best value
head(summary(ans8, order = value)) # best few
## head(summary(ans8, order = "value")) # best few -- alternative syntax

## order by value.  Within those values the same to 3 decimals order by fevals.
## summary(ans8, order = list(round(value, 3), fevals), par.select = FALSE)
summary(ans8, order = "list(round(value, 3), fevals)", par.select = FALSE)

## summary(ans8, order = rownames, par.select = FALSE) # order by method name
summary(ans8, order = "rownames", par.select = FALSE) # same

summary(ans8, order = NULL, par.select = FALSE) # use input order
## summary(ans8, par.select = FALSE) # same

ans8["value"]
