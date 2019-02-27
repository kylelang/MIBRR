### Title:    Test MIBRR Package
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2019-FEB-26

rm(list = ls(all = TRUE))

                                        #library(devtools)
                                        #install_github("kylelang/MIBRR/source/MIBRR", ref = "develop")
                                        #install_github("kylelang/SURF/source/SURF")

                                        #source("subroutines.R")

library(MIBRR)
library(SURF)
                                        #library(parallel)
                                        #library(HyperbolicDist)
                                        #library(LaplacesDemon)
                                        #library(monomvn)
                                        #library(rstan)

                                        #rstan_options(auto_write = TRUE)

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

## MCEM estimation:
mibenOut <- miben(data       = dat1,
                  iterations = c(30, 10),
                  targetVars = c("y", paste0("x", c(1 : 3))),
                  control    = list(savePpSams = TRUE)
                  )

