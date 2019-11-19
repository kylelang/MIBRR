### Title:    Test MIBRR Package
### Author:   Kyle M. Lang
### Created:  2018-FEB-06
### Modified: 2018-FEB-13

rm(list = ls(all = TRUE))

library(mitools)
library(psych)
library(MIBRR)
library(parallel)
library(devtools)

install_github("kylelang/SURF/source/SURF")
library(SURF)

###---------------------------------------------------------------------------###

### Prepare Data for Testing ###

data(bfi)
tmp <- na.omit(bfi)

ed.d           <- model.matrix(~factor(tmp$education))[ , -1]
colnames(ed.d) <-
    c("finish_hs", "some_college", "college_grad", "graduate_degree")

male            <- tmp$gender
male[male == 2] <- 0

cn   <- setdiff(colnames(bfi), c("gender", "education"))
bfi2 <- data.frame(tmp[ , cn], male, ed.d)

rownames(bfi2) <- NULL

targets  <- list(mar = paste0(c("A", "E"), rep(c(1 : 5), each = 2)),
                 mcar = NA,
                 mnar = NA)
pm       <- list(mar = 0.3)
snr      <- list(mar = 5)
marPreds <- c("age",
              "male",
              "finish_hs",
              "some_college",
              "college_grad",
              "graduate_degree")
cn       <- c(targets$mar, marPreds)

dat1 <- bfi2[sample(c(1 : nrow(bfi2)), 500), cn]
dat2 <- imposeMissData(data    = dat1,
                       targets = targets,
                       preds   = marPreds,
                       pm      = pm,
                       snr     = snr)$data 

###---------------------------------------------------------------------------###

### Check fully Bayesian estimation ###

mibenOut <- miben(data         = dat2,
                  targetVars   = targets$mar,
                  ignoreVars   = NULL,
                  sampleSizes  = c(1000, 1000),
                  doMcem       = FALSE,
                  lam1PriorPar = c(1.0, 0.1),
                  lam2PriorPar = c(1.0, 0.1),
                  verbose      = TRUE)

mibenImps <- getImpData(mibenOut, 100)

miblOut <- mibl(data         = dat2,
                targetVars   = targets$mar,
                ignoreVars   = NULL,
                sampleSizes  = c(1000, 1000),
                doMcem       = FALSE,
                lam1PriorPar = c(1.0, 0.1),
                verbose      = TRUE)

miblImps <- getImpData(miblOut, 100)

benOut <- ben(data         = bfi2,
              y            = "A1",
              X            = c(paste0("E", 1 : 5), paste0("O", 1 : 5)),
              sampleSizes  = c(500, 500),
              doMcem       = FALSE,
              lam1PriorPar = c(1.0, 0.1),
              lam2PriorPar = c(1.0, 0.1),
              verbose      = TRUE)

benPars <- getParams(benOut, "A1")
            
blOut <- bl(data         = bfi2,
            y            = "A1",
            X            = c(paste0("E", 1 : 5), paste0("O", 1 : 5)),
            sampleSizes  = c(500, 500),
            doMcem       = FALSE,
            lam1PriorPar = c(1.0, 0.1),
            verbose      = TRUE)

blPars <- getParams(blOut, "A1")

###---------------------------------------------------------------------------###

### Visualize Priors ###

nObs <- 500000

## Gamma:
shape <- 1.0
rate  <- 0.1

par(mfrow = c(4, 4))
vals <- c(0.01, 0.1, 1.0, 10.0)
for(i in vals)
    for(j in vals) {
        out1 <- rgamma(nObs, shape = i, rate = j)
        plot(density(out1), main = paste("Shape =", i, "Rate =", j))
    }
}

par(mfrow = c(1, 1))
out1 <- rgamma(nObs, shape = shape, rate = rate)
plot(density(out1))

## GIG:
par(mfrow = c(4, 4))
vals <- c(0.01, 0.1, 1.0, 10.0)
for(i in vals)
    for(j in vals) {
        out2 <- MIBRR:::drawGig(nObs, 1.0, i, j)
        plot(density(out2), main = paste("Chi =", i, "Psi =", j))
    }
}

par(mfrow = c(1, 1))
out2 <- MIBRR:::drawGig(nObs, 1.0, 1.0, 0.1)
plot(density(out2))

###---------------------------------------------------------------------------###

### Check Sensitivity to Priors ###

testFun1 <- function(rp, data, parms) {
    cat(paste0("Doing replication ", rp, "\n"))
    
    cn       <- parms$cn
    targets  <- parms$targets
    marPreds <- parms$marPreds
    pm       <- parms$pm
    snr      <- parms$snr
    nImps    <- parms$nImps
    keys     <- parms$keys
    
    dat1 <- data[sample(c(1 : nrow(data)), 500), cn]
    dat2 <- imposeMissData(data    = dat1,
                           targets = targets,
                           preds   = marPreds,
                           pm      = pm,
                           snr     = snr)$data

    outList <- list()
    for(v in 1 : nrow(parms$vals)) {
        vals <- parms$vals[v, ]
        if(parms$miben)
            fit <- miben(data         = dat2,
                         targetVars   = targets$mar,
                         ignoreVars   = NULL,
                         doMcem       = FALSE,
                         sampleSizes  = parms$samSizes,
                         lam1PriorPar = c(1.0, vals[1]),
                         lam2PriorPar = c(1.0, vals[2]),
                         verbose      = parms$verbose)

        else
            fit <- mibl(data         = dat2,
                        targetVars   = targets$mar,
                        ignoreVars   = NULL,
                        doMcem       = FALSE,
                        sampleSizes  = parms$samSizes,
                        lam1PriorPar = c(1.0, vals[1]),
                        verbose      = parms$verbose)
        imps <- getImpData(fit, nImps)
        
        resList <- list()
        for(m in 1 : nImps) {
            scores       <- scoreItems(keys  = keys, items = imps[[m]])$scores
            resList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                              mA = mean(scores[ , "agree"]),
                              mE = mean(scores[ , "extra"])
                              )
        }
        
        tag            <- ifelse(parms$miben,
                                 paste0("l1 = ", vals[1], " l2 = ", vals[2]),
                                 paste0("l1 = ", vals[1])
                                 )
        outList[[tag]] <- colMeans(do.call(rbind, resList))
    }
    outList
} # END testFun()


nReps <- 10
nImps <- 100
keys  <- list(agree = c("-A1", "A2", "A3", "A4", "A5"),
              extra = c("-E1", "-E2", "E3", "E4", "E5")
              )

parms          <- list()
parms$cn       <- cn
parms$targets  <- targets
parms$marPreds <- marPreds
parms$pm       <- pm
parms$snr      <- snr
parms$nImps    <- nImps
parms$keys     <- keys
parms$verbose  <- FALSE
parms$miben    <- TRUE
parms$samSizes <- c(500, 500)

if(parms$miben) {
    parms$vals <- expand.grid(c(0.01, 0.1, 0.5, 1.0, 5.0, 10.0),
                              c(0.01, 0.1, 0.5, 1.0, 5.0, 10.0)
                              )
} else {
    parms$vals <- matrix(c(0.01, 0.1, 0.5, 1.0, 5.0, 10.0))
}

simOut <- mclapply(X        = c(1 : nReps),
                   FUN      = testFun1,
                   data     = bfi2,
                   parms    = parms,
                   mc.cores = 4)

tmp  <- do.call(rbind, lapply(simOut, unlist))
tmp2 <- matrix(colMeans(tmp), ncol = 3, byrow = TRUE)
tmp3 <- apply(tmp2, 2, scale, scale = FALSE)

## How much Monte Carlo variability?
tmp4 <- matrix(apply(tmp, 2, sd), ncol = 3, byrow = TRUE)
tmp4

## How much prior-driven variability?
plot(tmp3[ , 1], ylim = range(tmp3), type = "l")
lines(tmp3[ , 2], col = "red")
lines(tmp3[ , 3], col = "blue")

###---------------------------------------------------------------------------###

### Compare Fully Bayesian and MCEM Estimation ###

testFun2 <- function(rp, data, parms) {
    cat(paste0("Doing replication ", rp, "\n"))
    
    cn       <- parms$cn
    targets  <- parms$targets
    marPreds <- parms$marPreds
    pm       <- parms$pm
    snr      <- parms$snr
    nImps    <- parms$nImps
    keys     <- parms$keys
    
    dat1 <- data[sample(c(1 : nrow(data)), 500), cn]
    dat2 <- imposeMissData(data    = dat1,
                           targets = targets,
                           preds   = marPreds,
                           pm      = pm,
                           snr     = snr)$data

    if(parms$miben) {
        mcemFit <- miben(data       = dat2,
                         targetVars = targets$mar,
                         ignoreVars = NULL,
                         iterations = parms$iterations,
                         verbose    = parms$verbose)
        mcemImps <- getImpData(mcemFit, nImps)
        
        bayesFit <- miben(data         = dat2,
                          targetVars   = targets$mar,
                          ignoreVars   = NULL,
                          doMcem       = FALSE,
                          sampleSizes  = parms$samSizes,
                          lam1PriorPar = parms$l1Par,
                          lam2PriorPar = parms$l2Par,
                          verbose      = parms$verbose)
        bayesImps <- getImpData(bayesFit, nImps)
    }
    else {
        mcemFit <- mibl(data       = dat2,
                        targetVars = targets$mar,
                        ignoreVars = NULL,
                        iterations = parms$iterations,
                        verbose    = parms$verbose)
        mcemImps <- getImpData(mcemFit, nImps)
        
        bayesFit <- mibl(data         = dat2,
                         targetVars   = targets$mar,
                         ignoreVars   = NULL,
                         doMcem       = FALSE,
                         sampleSizes  = parms$samSizes,
                         lam1PriorPar = parms$l1Par,
                         verbose      = parms$verbose)
        bayesImps <- getImpData(bayesFit, nImps)
    }
    
    mcemRes <- bayesRes <- list()
    for(m in 1 : nImps) {
        ## MCEM estimates:
        scores       <- scoreItems(keys  = keys, items = mcemImps[[m]])$scores
        mcemRes[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                          mA = mean(scores[ , "agree"]),
                          mE = mean(scores[ , "extra"])
                          )
        ## Bayes estimates:
        scores        <- scoreItems(keys  = keys, items = bayesImps[[m]])$scores
        bayesRes[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                           mA = mean(scores[ , "agree"]),
                           mE = mean(scores[ , "extra"])
                           )
    }
    
    list(mcem  = colMeans(do.call(rbind, mcemRes)),
         bayes = colMeans(do.call(rbind, bayesRes))
         )
} # END testFun()


nReps <- 8
nImps <- 100
keys  <- list(agree = c("-A1", "A2", "A3", "A4", "A5"),
              extra = c("-E1", "-E2", "E3", "E4", "E5")
              )

parms            <- list()
parms$cn         <- cn
parms$targets    <- targets
parms$marPreds   <- marPreds
parms$pm         <- pm
parms$snr        <- snr
parms$nImps      <- nImps
parms$keys       <- keys
parms$verbose    <- FALSE
parms$miben      <- TRUE
parms$samSizes   <- c(500, 500)
parms$iterations <- c(50, 10)
parms$l1Par      <- c(1.0, 0.1)
parms$l2Par      <- c(1.0, 0.1)

simOut <- mclapply(X        = c(1 : nReps),
                   FUN      = testFun2,
                   data     = bfi2,
                   parms    = parms,
                   mc.cores = 4)

tmp  <- do.call(rbind, lapply(simOut, unlist))

tmp

mcemFrame  <- tmp[ , grep("mcem", colnames(tmp))]
bayesFrame <- tmp[ , grep("bayes", colnames(tmp))]

## Complete data result: 
scores  <- scoreItems(keys = keys, items = bfi2)$scores
compRes <- c(r  = cor(scores[ , 1], scores[ , 2]),
             mA = mean(scores[ , "agree"]),
             mE = mean(scores[ , "extra"])
             )

## Difference in two methods (as proportion of complete-data result):
100 * colMeans(mcemFrame - bayesFrame) / compRes

## Percent Relative Bias:
100 * (colMeans(mcemFrame) - compRes) / compRes
100 * (colMeans(bayesFrame) - compRes) / compRes

# Monte Carlo SD:
apply(mcemFrame, 2, sd)
apply(bayesFrame, 2, sd)
