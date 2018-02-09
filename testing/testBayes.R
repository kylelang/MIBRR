### Title:    Test MIBRR Package
### Author:   Kyle M. Lang
### Created:  2018-FEB-06
### Modified: 2018-FEB-09

rm(list = ls(all = TRUE))

                                        #install.packages("HyperbolicDist", repos = "http://cloud.r-project.org")

library(mitools)
library(psych)
library(MIBRR)
library(devtools)
library(parallel)
library(MCMCpack)
library(statmod)
library(HyperbolicDist)

install_github("kylelang/SURF/source/SURF")
library(SURF)

saveDate <- format(Sys.time(), "%Y%m%d")

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

### Check fully Bayesian estimation of MIBEN and MIBL ###

mibenOut <- miben(data         = dat2,
                  targetVars   = targets$mar,
                  ignoreVars   = NULL,
                  sampleSizes  = c(1000, 1000),
                  doMcem       = FALSE,
                  lam1PriorPar = c(1.0, 0.1),
                  lam2PriorPar = c(1.0, 0.1),
                  verbose      = TRUE)

mibenImps <- MIBRR::complete(mibenOut, 100)


miblOut <- mibl(data         = dat2,
                targetVars   = targets$mar,
                ignoreVars   = NULL,
                sampleSizes  = c(1000, 1000),
                doMcem       = FALSE,
                lam1PriorPar = c(1.0, 0.1),
                verbose      = TRUE)

miblImps <- MIBRR::complete(miblOut, 100)


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

testFun <- function(rp, data, parms) {
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
        imps <- MIBRR::complete(fit, nImps)
        
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


nReps <- 20
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
parms$miben    <- FALSE
parms$samSizes <- c(500, 500)
parms$vals     <- matrix(c(0.01, 0.1, 0.5, 1.0, 5.0, 10.0))

                                        #parms$vals     <- expand.grid(c(0.01, 0.1, 0.5, 1.0, 5.0, 10.0),
                                        #                              c(0.01, 0.1, 0.5, 1.0, 5.0, 10.0)
                                        #                              )

simOut <- mclapply(X        = c(1 : nReps),
                   FUN      = testFun,
                   data     = bfi2,
                   parms    = parms,
                   mc.cores = 4)

tmp  <- do.call(rbind, lapply(simOut, unlist))
tmp2 <- matrix(colMeans(tmp), ncol = 3, byrow = TRUE)
tmp3 <- apply(tmp2, 2, scale, scale = FALSE)

tmp4 <- matrix(apply(tmp, 2, sd), ncol = 3, byrow = TRUE)

plot(tmp3[ , 1], ylim = range(tmp3), type = "l")
lines(tmp3[ , 2], col = "red")
lines(tmp3[ , 3], col = "blue")

mibenFrame <- tmp[ , grep("miben", colnames(tmp))]
miblFrame  <- tmp[ , grep("mibl", colnames(tmp))]
miceFrame  <- tmp[ , grep("mice", colnames(tmp))]

## Complete data result: 
scores  <- scoreItems(keys = keys, items = bfi2)$scores
compRes <- c(r  = cor(scores[ , 1], scores[ , 2]),
             mA = mean(scores[ , "agree"]),
             mE = mean(scores[ , "extra"])
             )

## Percent Relative Bias:
100 * (colMeans(mibenFrame) - compRes) / compRes
100 * (colMeans(miblFrame) - compRes) / compRes
100 * (colMeans(miceFrame) - compRes) / compRes

# Monte Carlo SD:
apply(mibenFrame, 2, sd)
apply(miblFrame, 2, sd)
apply(miceFrame, 2, sd)
