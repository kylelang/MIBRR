### Title:    MIBRR Imputation Simulation
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2018-NOV-08

rm(list = ls(all = TRUE))

library(mitools)
library(psych)
library(MIBRR)
library(parallel)
library(SURF)

###--------------------------------------------------------------------------###

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

###--------------------------------------------------------------------------###

### Define a function to run each replication ###

testFun <- function(rp, data, env) {
    cn       <- env$cn
    targets  <- env$targets
    marPreds <- env$marPreds
    pm       <- env$pm
    snr      <- env$snr
    nImps    <- env$nImps
    keys     <- env$keys
    
    dat1 <- data[sample(c(1 : nrow(data)), 500), cn]
    dat2 <- imposeMissData(data    = dat1,
                           targets = targets,
                           preds   = marPreds,
                           pm      = pm,
                           snr     = snr)$data

    mibenOut <- miben(data       = dat2,
                      targetVars = targets$mar,
                      ignoreVars = NULL,
                      iterations = c(50, 10),
                      verbose    = FALSE)
    mibenImps <- getImpData(mibenOut, nImps)
    
    miblOut <- mibl(data       = dat2,
                    targetVars = targets$mar,
                    ignoreVars = NULL,
                    iterations = c(50, 10),
                    verbose    = FALSE)
    miblImps <- getImpData(miblOut, nImps)

    vanOut <- vanilla(data       = dat2,
                      targetVars = targets$mar,
                      ignoreVars = NULL,
                      verbose    = FALSE)
    vanImps <- getImpData(vanOut, nImps)
    
    miceOut <-
        mice(dat2, m = nImps, maxit = 10, method = "norm", printFlag = FALSE)

    mibenList <- miblList <- vanList <- miceList <- list()
    for(m in 1 : nImps) {
        ## MIBEN estimates:
        scores         <- scoreItems(keys  = keys,
                                     items = mibenImps[[m]])$scores
        mibenList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                            mA = mean(scores[ , "agree"]),
                            mE = mean(scores[ , "extra"])
                            )

        ## MIBL estimates:
        scores        <- scoreItems(keys  = keys,
                                    items = miblImps[[m]])$scores
        miblList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                           mA = mean(scores[ , "agree"]),
                           mE = mean(scores[ , "extra"])
                           )

        ## Vanilla estimates:
        scores         <- scoreItems(keys  = keys,
                                     items = vanImps[[m]])$scores
        vanList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                          mA = mean(scores[ , "agree"]),
                          mE = mean(scores[ , "extra"])
                          )

        ## MICE estimates:
        scores        <- scoreItems(keys  = keys,
                                    items = mice::complete(miceOut, m))$scores
        miceList[[m]] <- c(r  = cor(scores[ , 1], scores[ , 2]),
                           mA = mean(scores[ , "agree"]),
                           mE = mean(scores[ , "extra"])
                           )
    }

    list(miben = colMeans(do.call(rbind, mibenList)),
         mibl  = colMeans(do.call(rbind, miblList)),
         van   = colMeans(do.call(rbind, vanList)),
         mice  = colMeans(do.call(rbind, miceList))
         )
} # END testFun()

###--------------------------------------------------------------------------###

### Run the simulation ###

nReps <- 100
nImps <- 20
keys  <- list(agree = c("-A1", "A2", "A3", "A4", "A5"),
              extra = c("-E1", "-E2", "E3", "E4", "E5")
              )

simOut <- mclapply(X        = c(1 : nReps),
                   FUN      = testFun,
                   data     = bfi2,
                   env      = parent.frame(),
                   mc.cores = 4)

###--------------------------------------------------------------------------###

### Pool the results ###

tmp <- do.call(rbind, lapply(simOut, unlist))

mibenFrame <- tmp[ , grep("miben", colnames(tmp))]
miblFrame  <- tmp[ , grep("mibl", colnames(tmp))]
vanFrame   <- tmp[ , grep("van", colnames(tmp))]
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
100 * (colMeans(vanFrame) - compRes) / compRes
100 * (colMeans(miceFrame) - compRes) / compRes
