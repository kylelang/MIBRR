### Title:    Build R Packages
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2016-MAY-05
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))
require(RcppEigen)

system("rm -r packageSource/mibrr/src/nlopt/* \
        rm packageSource/mibrr/src/RcppExports.cpp \
        rm packageSource/mibrr/R/RcppExports.R \
        rm packageSource/mibrr/src/*.o packageSource/mibrr/src/*.so")
system("cp -r ~/data/software/miscPackages/nlopt-2.4.2/* packageSource/mibrr/src/nlopt/")
Rcpp::compileAttributes("packageSource/mibrr")
install.packages("packageSource/mibrr", repos = NULL, type = "source")

library(mibrr)

data(mibrrExampleData)

testOut <- miben(rawData      = mibrrExampleData,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                 returnParams = TRUE)

testOut2 <- mibl(rawData      = mibrrExampleData,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                 returnParams = TRUE)


?read.dcf

mibrrL()
mibrrW()
