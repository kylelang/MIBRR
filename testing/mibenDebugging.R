### Title:    Debug MIBRR Package
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2019-JAN-15

rm(list = ls(all = TRUE))

library(MIBRR)

data(mibrrExampleData)

undebug(miben)

## MCEM estimation:
mibenOut <- miben(data       = mibrrExampleData,
                  iterations = c(30, 10),
                  targetVars = c("y", paste0("x", c(1 : 3))),
                  ignoreVars = "idNum",
                  control = list(optTraceLevel = 1,
                                 optMethod = c("Nelder-Mead", "BFGS"),
                                 optBoundLambda = FALSE
                                 )
                  )

?miben

## Fully Bayesian estimation:
mibenOut <- miben(data         = mibrrExampleData,
                  targetVars   = c("y", paste0("x", c(1 : 3))),
                  ignoreVars   = "idNum",
                  sampleSizes  = c(500, 500),
                  doMcem       = FALSE,
                  lam1PriorPar = c(1.0, 0.1),
                  lam2PriorPar = c(1.0, 0.1)
                  )

?optimx
?mibl

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
