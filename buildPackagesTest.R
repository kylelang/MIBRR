### Title:    Build R Packages
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2016-MAY-05
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))

install.packages(c("mitools"),
                 repos = "http://rweb.quant.ku.edu/cran")

library(RcppEigen)

system("rm -r packageSource/mibrr/src/nlopt/* \
        rm packageSource/mibrr/src/RcppExports.cpp \
        rm packageSource/mibrr/R/RcppExports.R \
        rm packageSource/mibrr/src/*.o packageSource/mibrr/src/*.so")
system("cp -r ~/data/software/miscPackages/nlopt-2.4.2/* packageSource/mibrr/src/nlopt/")
Rcpp::compileAttributes("packageSource/mibrr")
install.packages("packageSource/mibrr", repos = NULL, type = "source")

library(mitools)
library(mibrr)

data(mibrrExampleData)

debug(miben)
undebug(miben)

testOut <- miben(rawData      = mibrrExampleData,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                 returnParams = TRUE)

fitOut <- lapply(testOut$imps,
                 FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                 )
MIcombine(fitOut)

testOut2 <- mibl(rawData      = mibrrExampleData,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                 returnParams = TRUE)

fitOut2 <- lapply(testOut$imps,
                 FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                 )
MIcombine(fitOut2)
