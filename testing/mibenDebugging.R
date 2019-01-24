### Title:    Debug MIBRR Package
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2019-JAN-21

rm(list = ls(all = TRUE))

library(MIBRR)

data(mibrrExampleData)

## MCEM estimation:
mibenOut <- miben(data       = mibrrExampleData,
                  iterations = c(30, 10),
                  targetVars = c("y", paste0("x", c(1 : 3))),
                  ignoreVars = "idNum",
                  control = list(optTraceLevel = 0,
                                 optMethod = "L-BFGS-B", #c("Nelder-Mead", "BFGS"),
                                 optBoundLambda = TRUE
                                 )
                  )

getParams(mibenOut, "y")

## Fully Bayesian estimation:
mibenOut <- miben(data         = mibrrExampleData,
                  targetVars   = c("y", paste0("x", c(1 : 3))),
                  ignoreVars   = "idNum",
                  sampleSizes  = c(500, 500),
                  doMcem       = FALSE,
                  lam1PriorPar = c(1.0, 0.1),
                  lam2PriorPar = c(1.0, 0.1)
                  )

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


missList   <- readRDS("missList.rds")
missCounts <- readRDS("missCounts.rds")
dat1       <- readRDS("data.rds")

respCounts <- nrow(dat1) - missCounts

X1 <- MIBRR:::getX(data        = as.matrix(dat1),
                   missIndices = missList,
                   respCounts  = respCounts,
                   noMiss      = FALSE,
                   xOnly       = TRUE,
                   obsY        = TRUE,
                   scale       = TRUE,
                   targetIndex = 0)

X2 <- MIBRR:::getX(data        = as.matrix(dat1),
                   missIndices = missList,
                   respCounts  = respCounts,
                   noMiss      = FALSE,
                   xOnly       = TRUE,
                   obsY        = FALSE,
                   scale       = TRUE,
                   targetIndex = 0)

X3 <- MIBRR:::getX(data        = as.matrix(dat1),
                   missIndices = missList,
                   respCounts  = respCounts,
                   noMiss      = FALSE,
                   xOnly       = TRUE,
                   obsY        = TRUE,
                   scale       = FALSE,
                   targetIndex = 0)

y <- MIBRR:::getY(data        = as.matrix(dat1),
                  missIndices = missList,
                  respCounts  = respCounts,
                  noMiss      = FALSE,
                  targetIndex = 0)

y
mean(y)
var(y)
sd(y)
crossprod(y)

X1
colMeans(X1)
apply(X1, 2, sd)
apply(X1, 2, var)
apply(X1, 2, crossprod)

X2
colMeans(X2)
apply(X2, 2, sd)
apply(X2, 2, var)
apply(X2, 2, crossprod)

