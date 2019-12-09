### Title:    Test MIBRR Package
### Author:   Kyle M. Lang
### Created:  2014-12-07
### Modified: 2019-12-09

rm(list = ls(all = TRUE))

                                        #library(devtools)
                                        #install_github("kylelang/MIBRR/source/MIBRR", ref = "multithread")
                                        #install_github("kylelang/SURF/source/SURF")

source("testingSubroutines.R")

library(MIBRR)
library(SURF)
library(parallel)
library(HyperbolicDist)
library(LaplacesDemon)
library(monomvn)
library(rstan)

rstan_options(auto_write = TRUE)

###--------------------------------------------------------------------------###

## Generate some data:
dat0 <- simRegData(nObs  = 500,
                   nVars = 10,
                   r2    = 0.5,
                   sigma = 0.2,
                   beta  = matrix(c(0.25, rep(0.75, 10)))
                   )

## Impose missing values:
dat1 <- imposeMissData(data    = dat0,
                       targets = list(mar  = colnames(dat0)[1 : 8],
                                      mcar = colnames(dat0)[9 : 11]
                                      ),
                       preds   = colnames(dat0)[9 : 11],
                       pm      = list(mar = 0.3, mcar = 0.1),
                       snr     = 5.0)$data

###--------------------------------------------------------------------------###

### Run Internal (Unit-ish) Tests ###

MIBRR:::testMissIndex(dat1)
MIBRR:::testDataProcessing(dat1)
MIBRR:::testDataScaling(dat1)
MIBRR:::testMissFill(dat1)
MIBRR:::testSamplers()

###--------------------------------------------------------------------------###

### Compare MIBRR::bl to a reference implementation ###

xNames <- setdiff(colnames(dat0), "y")
iters  <- c(500, 20)
sams   <- list(c(25, 25), c(100, 100), c(1000, 1000))

### MCEM Estimation:

## Estimate the model using monomvn::blasso:
bl0 <- bl0Mcem(data      = dat0,
               yName     = "y",
               xNames    = xNames,
               iters     = iters,
               sams      = sams,
               intercept = TRUE,
               norm      = TRUE)

## Estimate the model using MIBRR:bl:
bl1 <-
    bl(data        = dat0,
       y           = "y",
       X           = xNames,
       iterations  = iters,
       sampleSizes = sams,
       verbose     = TRUE,
       control     = list(lambda1Starts = 1.0, useBetaMeans = FALSE)
       )

## Generate posterior predictive samples:
pp0 <- predBl0(bl0$out, X = as.matrix(dat0[ , xNames]))[ , -1]
pp1 <- postPredict(bl1, newData = dat0, nDraws = sams[[3]][2])[[1]]

## Plot 25 randomly sample posterior predictive densities:
par(mfrow = c(5, 5))

for(x in sample(1 : nrow(dat0), 25)) {
    d0 <- density(pp0[x, ])
    d1 <- density(pp1[x, ])

    plot(d0, ylim = range(d0$y, d1$y), xlim = range(d0$x, d1$x))
    lines(d1, col = "red")
}

## Extract parameter samples:
b0 <- cbind(bl0$out$mu[-1], bl0$out$beta[-1, ])
b1 <- getParams(bl1, "y")$beta

par(mfrow = c(4, 4))

for(x in 1 : ncol(b1)) {
    d0 <- density(b0[ , x])
    d1 <- density(b1[ , x])
    
    plot(d0, ylim = range(d0$y, d1$y), xlim = range(d0$x, d1$x))
    lines(d1, col = "red")
}   

### Fully Bayesian Estimation:

## Estimate the model using monomvn::blasso:
bl0 <- blasso(y         = dat0$y,
              X         = dat0[ , xNames],
              T         = sams[[3]][2],
              thin      = sams[[3]][1],
              lambda2   = 1.0,
              s2        = with(dat0, var(y - mean(y))),
              beta      = rnorm(ncol(dat0) - 1),
              rd        = c(0.5, 0.5),
              RJ        = FALSE,
              rao.s2    = FALSE,
              icept     = TRUE,
              normalize = TRUE)

## Estimate the model using MIBRR:bl:
bl1 <-
    bl(data         = dat0,
       y            = "y",
       X            = xNames,
       sampleSizes  = sams[[3]],
       verbose      = TRUE,
       doMcem       = FALSE,
       lam1PriorPar = c(0.5, 0.5),
       control      = list(lambda1Starts = 1.0, useBetaMeans = FALSE)
       )

## Generate posterior predictive samples:
pp0 <- predBl0(bl0, X = as.matrix(dat0[ , xNames]))[ , -1]
pp1 <-
    postPredict(bl1, newData = dat0, nDraws = sams[[3]][2], scale = FALSE)[[1]]

## Plot 25 randomly sample posterior predictive densities:
par(mfrow = c(5, 5))

for(x in sample(1 : nrow(dat0), 25)) {
    d0 <- density(pp0[x, ])
    d1 <- density(pp1[x, ])

    plot(d0, ylim = range(d0$y, d1$y), xlim = range(d0$x, d1$x))
    lines(d1, col = "red")
}

## Extract and visualize sigma samples:
s0 <- bl0$s2
s1 <- getParams(bl1, "y")$sigma

par(mfrow = c(1, 1))

d0 <- density(s0)
d1 <- density(s1)

plot(d0, ylim = range(d0$y, d1$y), xlim = range(d0$x, d1$x))
lines(d1, col = "red")

## Extract parameter samples:

b0 <- cbind(bl0$mu[-1], bl0$beta[-1, ])
b1 <- getParams(bl1, "y")$beta

par(mfrow = c(4, 4))

for(x in 1 : ncol(b1)) {
    d0 <- density(b0[ , x])
    d1 <- density(b1[ , x])
    
    plot(d0, ylim = range(d0$y, d1$y), xlim = range(d0$x, d1$x))
    lines(d1, col = "red")
}   

###--------------------------------------------------------------------------###

### Compare MIBRR::ben to a reference implementation ###

## Compile the Stan model:
benMod <- stan_model("ben_model.stan")

## Create a data list for Stan:
stanData <- list(N = nrow(dat0),
                 P = ncol(dat0) - 1,
                 y = dat0$y,
                 X = dat0[ , xNames])

## Sample from the Stan model:
ben0 <- sampling(object = benMod,
                 data   = stanData,
                 iter   = sum(sams[[3]]),
                 warmup = sams[[3]][1],
                 chains = 1)

## Estimate the model using MIBRR:ben:
ben1 <-
    ben(data         = dat0,
        y            = "y",
        X            = xNames,
        iterations   = iters,
        sampleSizes  = sams[[3]],
        doMcem       = FALSE,
        lam1PriorPar = c(1.0, 0.1),
        lam2PriorPar = c(1.0, 0.1),
        verbose      = TRUE,
        control      = list(useBetaMeans = FALSE)
        )

## Generate posterior predictive samples:
pp0 <- predBen0(obj = ben0, X = as.matrix(dat0[ , xNames]))
pp1 <-
    postPredict(ben1, newData = dat0, nDraws = sams[[3]][2], scale = FALSE)[[1]]

## Plot 25 randomly sample posterior predictive densities:
par(mfrow = c(5, 5))

for(x in sample(1 : nrow(dat0), 25)) {
    d0 <- density(pp0[x, ])
    d1 <- density(pp1[x, ])
    
    plot(d0, ylim = range(d0$y, d1$y), xlim = range(d0$x, d1$x))
    lines(d1, col = "red")
}

## Extract and visualize sigma samples:
s0 <- extract(ben0, pars = "sigma")[[1]]
s1 <- getParams(ben1, "y")$sigma

par(mfrow = c(1, 1))

d0 <- density(s0^2)
d1 <- density(s1)

plot(d0, ylim = range(d0$y, d1$y), xlim = range(d0$x, d1$x))
lines(d1, col = "red")

## Extract parameter samples:
tmp <- extract(ben0, pars = c("mu", "beta"))
b0  <- cbind(tmp[["mu"]], tmp[["beta"]])
b1  <- getParams(ben1, "y")$beta

par(mfrow = c(4, 4))

for(x in 1 : ncol(b1)) {
    d0 <- density(b0[ , x])
    d1 <- density(b1[ , x])
    
    plot(d0, ylim = range(d0$y, d1$y), xlim = range(d0$x, d1$x))
    lines(d1, col = "red")
}   

###--------------------------------------------------------------------------###

### Compare MIBRR::bvr to a reference implementation ###

## Run both versions:
par0 <- bReg(data = dat0, y = "y", X = xNames, nSams = 10000)
br1  <- bvr(data        = dat0,
            y           = "y",
            X           = xNames,
            sampleSizes = rep(10000, 2),
            ridge       = 0.0)
par1 <- getParams(br1, "y")

## Extract and visualize sigma samples:
s0 <- par0$sigma2
s1 <- par1$sigma

par(mfrow = c(1, 1))

d0 <- density(s0)
d1 <- density(s1)

plot(d0, ylim = range(d0$y, d1$y), xlim = range(d0$x, d1$x))
lines(d1, col = "red")
   
## Extract and visualize beta samples:
b0 <- par0$beta
b1 <- getParams(br1, "y")$beta

par(mfrow = c(3, 4))

for(x in 1 : ncol(b1)) {
    d0 <- density(b0[ , x])
    d1 <- density(b1[ , x])
    
    plot(d0, ylim = range(d0$y, d1$y), xlim = range(d0$x, d1$x))
    lines(d1, col = "red")
}   

###--------------------------------------------------------------------------###

### Do Posterior Predictive Checks ###

## BEN ##

## MCEM estimation:
benOut <- ben(data       = dat0,
              y          = "y",
              X          = paste0("x", c(1 : 3)),
              iterations = c(30, 10),
              control    = list(savePpSams = TRUE)
              )

par(mfrow = c(1, 1))
MIBRR:::ppCheck(benOut)

## Fully Bayesian estimation:
benOut <- ben(data         = dat0,
              y            =  "y",
              X            = paste0("x", c(1 : 3)),
              sampleSizes  = c(500, 500),
              doMcem       = FALSE,
              lam1PriorPar = c(1.0, 0.1),
              lam2PriorPar = c(1.0, 0.1),
              control      = list(savePpSams = TRUE)
              )

par(mfrow = c(1, 1))
MIBRR:::ppCheck(benOut)

## MIBEN ##

## MCEM estimation:
mibenOut <- miben(data       = dat1,
                  iterations = c(30, 10),
                  targetVars = c("y", paste0("x", c(1 : 3))),
                  control    = list(savePpSams = TRUE)
                  )

par(mfrow = c(2, 2))
MIBRR:::ppCheck(mibenOut)

## Fully Bayesian estimation:
mibenOut <- miben(data         = dat1,
                  targetVars   = c("y", paste0("x", c(1 : 3))),
                  sampleSizes  = c(500, 500),
                  doMcem       = FALSE,
                  lam1PriorPar = c(1.0, 0.1),
                  lam2PriorPar = c(1.0, 0.1),
                  control      = list(savePpSams = TRUE)
                  )

par(mfrow = c(2, 2))
MIBRR:::ppCheck(mibenOut)

## BL ##

## MCEM estimation:
blOut <- bl(data       = dat0,
            y          = "y",
            X          = paste0("x", c(1 : 3)),
            iterations = c(30, 10),
            control    = list(savePpSams = TRUE)
            )

par(mfrow = c(1, 1))
MIBRR:::ppCheck(blOut)

## Fully Bayesian estimation:
blOut <- bl(data         = dat0,
            y            =  "y",
            X            = paste0("x", c(1 : 3)),
            sampleSizes  = c(500, 500),
            doMcem       = FALSE,
            lam1PriorPar = c(1.0, 0.1),
            control      = list(savePpSams = TRUE)
            )

par(mfrow = c(1, 1))
MIBRR:::ppCheck(blOut)

## MIBL ##

## MCEM estimation:
miblOut <- mibl(data       = dat1,
                iterations = c(30, 10),
                targetVars = c("y", paste0("x", c(1 : 3))),
                control    = list(savePpSams = TRUE)
                )

par(mfrow = c(2, 2))
MIBRR:::ppCheck(miblOut)

## Fully Bayesian estimation:
miblOut <- mibl(data         = dat1,
                targetVars   = c("y", paste0("x", c(1 : 3))),
                sampleSizes  = c(500, 500),
                doMcem       = FALSE,
                lam1PriorPar = c(1.0, 0.1),
                control      = list(savePpSams = TRUE)
                )

par(mfrow = c(2, 2))
MIBRR:::ppCheck(miblOut)

## Vanilla Bayesian Regression ##

bvrOut <- bvr(data        = dat0,
              y           = "y",
              X           = paste0("x", c(1 : 3)),
              sampleSizes = c(500, 500),
              control     = list(savePpSams = TRUE)
              )

par(mfrow = c(1, 1))
MIBRR:::ppCheck(bvrOut)

## Vanilla MI ##

vanOut <- vanilla(data        = dat1,
                  targetVars  = c("y", paste0("x", c(1 : 3))),
                  sampleSizes = c(500, 500),
                  ridge       = 0.0,
                  control     = list(savePpSams = TRUE)
                  )

par(mfrow = c(2, 2))
MIBRR:::ppCheck(vanOut)

###--------------------------------------------------------------------------###

### Graphically Test Random Variate Samplers ###

## MVN sampler:
nObs           <- 500000
mvnMu          <- rep(10, 3)
mvnSigma       <- matrix(5, 3, 3)
diag(mvnSigma) <- 20

out1.1 <- MIBRR:::drawMvn(nObs, mvnMu, mvnSigma, 235711)
out1.2 <- rmvnorm(nObs, mvnMu, mvnSigma)

par(mfrow = c(1, 3))
for(i in 1 : length(mvnMu)) {
    plot(density(out1.1[ , i]), col = "red")
    lines(density(out1.2[ , i]), col = "blue")
}

## Inverse gamma sampler:
gamShape <- 10
gamScale <- 10

out2.1 <- MIBRR:::drawInvGamma(nObs, gamShape, gamScale, 235711)
out2.2 <- rinvgamma(nObs, gamShape, gamScale)

plot(density(out2.1), col = "red")
lines(density(out2.2), col = "blue")

## Inverse gaussian sampler:
igMu  <- 1
igLam <- 2

out3.1 <- MIBRR:::drawInvGauss(nObs, igMu, igLam, 235711)
out3.2 <- rinvgaussian(nObs, igMu, igLam)

plot(density(out3.1), col = "red")
lines(density(out3.2), col = "blue")

## GIG sampler:
gigLam <- 1
gigChi <- 2
gigPsi <- 2

out4.1 <- MIBRR:::drawGig(nObs, gigLam, gigChi, gigPsi, 235711)
out4.2 <- rgig(nObs, c(gigLam, gigChi, gigPsi))

plot(density(out4.1), col = "red")
lines(density(out4.2), col = "blue")

## Scaled inverse chi-squared sampler:
df    <- 100
scale <- 10

out5.1 <- MIBRR:::drawInvChiSq(nObs, df, scale, 235711)
out5.2 <- rinvchisq(nObs, df, scale)

plot(density(out5.1), col = "red")
lines(density(out5.2), col = "blue")

###--------------------------------------------------------------------------###

### Check Documentation Examples ###

## MIBEN ##
?miben

data(mibrrExampleData)

## MCEM estimation:
mibenOut <- miben(data       = mibrrExampleData,
                  iterations = c(30, 10),
                  targetVars = c("y", paste0("x", c(1 : 3))),
                  ignoreVars = "idNum")

## Fully Bayesian estimation:
mibenOut <- miben(data         = mibrrExampleData,
                  targetVars   = c("y", paste0("x", c(1 : 3))),
                  ignoreVars   = "idNum",
                  sampleSizes  = c(500, 500),
                  doMcem       = FALSE,
                  lam1PriorPar = c(1.0, 0.1),
                  lam2PriorPar = c(1.0, 0.1)
                  )

## MIBL ##
?mibl

data(mibrrExampleData)

## MCEM estimation:
miblOut <- mibl(data       = mibrrExampleData,
                iterations = c(50, 10),
                targetVars = c("y", paste0("x", c(1 : 3))),
                ignoreVars = "idNum")

## Fully Bayesian estimation:
miblOut <- mibl(data         = mibrrExampleData,
                targetVars   = c("y", paste0("x", c(1 : 3))),
                ignoreVars   = "idNum",
                sampleSizes  = c(500, 500),
                doMcem       = FALSE,
                lam1PriorPar = c(1.0, 0.1)
                )

## BEN ##
?ben

data(mibrrExampleData)

## MCEM estimation:
benOut <- ben(data       = mibrrExampleData,
              y          = "y",
              X          = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
              iterations = c(30, 10)
              )

## Fully Bayesian estimation:
benOut <- ben(data         = mibrrExampleData,
              y            = "y",
              X            = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
              doMcem       = FALSE,
              sampleSizes  = c(500, 500),
              lam1PriorPar = c(1.0, 0.1),
              lam2PriorPar = c(1.0, 0.1)
              )
     
## BL ##
?bl

data(mibrrExampleData)

## MCEM estimation:
blOut <- bl(data       = mibrrExampleData,
            y          = "y",
            X          = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
            iterations = c(50, 10)
            )

## Fully Bayesian estimation:
blOut <- bl(data         = mibrrExampleData,
            y            = "y",
            X            = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
            doMcem       = FALSE,
            sampleSizes  = c(500, 500),
            lam1PriorPar = c(1.0, 0.1)
            )
     
## postPredict ##
?postPredict

data(predictData)

## Fit a Bayesian elastic net model:
benOut <- ben(data       = predictData$train,
              y          = "agree",
              X          = setdiff(colnames(predictData$train), "agree"),
              iterations = c(30, 10)
              )

## Generate 10 posterior predictions for 'y':
benPred <- postPredict(mibrrFit = benOut,
                       newData  = predictData$test,
                       nDraws   = 10)

## Generate posterior predictions for 'y' using MAP scores (the default):
benPred <- postPredict(mibrrFit = benOut,
                       newData  = predictData$test,
                       nDraws   = 0)

## Generate posterior predictions for 'y' using EAP scores:
benPred <- postPredict(mibrrFit = benOut,
                       newData  = predictData$test,
                       nDraws   = -1)

## Fit chained Bayesian elastic net models to support MI:
mibenOut <- miben(data = predictData$incomplete, iterations = c(30, 10))

## Generate posterior predictions from elementary imputation models:
mibenPred <- postPredict(mibrrFit = mibenOut, newData = predictData$test)

## getImpData ##
?getImpData

data(mibrrExampleData)

## Fit a Bayesian elastic net model:
mibenOut <- miben(data         = mibrrExampleData,
                  targetVars   = c("y", paste0("x", c(1 : 3))),
                  ignoreVars   = "idNum",
                  sampleSizes  = c(500, 500),
                  doMcem       = FALSE,
                  lam1PriorPar = c(1.0, 0.1),
                  lam2PriorPar = c(1.0, 0.1)
                  )

## Generate 25 imputed datasets:
mibenImps <- getImpData(mibrrFit = mibenOut, nImps = 25)

## getParams ##
?getParams

data(mibrrExampleData)

## Fit a Bayesian elastic net model:
benOut <- ben(data         = mibrrExampleData,
              y            = "y",
              X            = paste0("x", c(1 : 3)),
              sampleSizes  = c(500, 500),
              doMcem       = FALSE,
              lam1PriorPar = c(1.0, 0.1),
              lam2PriorPar = c(1.0, 0.1)
              )

## Extract posterior parameter samples:
benPars <- getParams(mibrrFit = benOut, target = "y")

## getField ##
?getField

data(mibrrExampleData)

## Fit a Bayesian elastic net model:
benOut <- ben(data         = mibrrExampleData,
              y            = "y",
              X            = paste0("x", c(1 : 3)),
              sampleSizes  = c(500, 500),
              doMcem       = FALSE,
              lam1PriorPar = c(1.0, 0.1),
              lam2PriorPar = c(1.0, 0.1)
              )

## Extract the potential scale reduction factors:
benRHats <- getField(mibrrFit = benOut, what = "rHats")

## Info Functions ##
?mibrrL
mibrrL()

?mibrrW
mibrrW()
