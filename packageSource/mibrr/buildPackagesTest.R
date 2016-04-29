### Title:    Build R Packages
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2016-APR-28
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))
require(RcppEigen)

system("rm mibrr/src/RcppExports.cpp mibrr/R/RcppExports.R mibrr/src/*.o mibrr/src/*.so")
Rcpp::compileAttributes("mibrr")
install.packages("mibrr", repos = NULL, type = "source")

library(mibrr)

data(mibrrExampleData)

testOut <- miben(rawData = mibrrExampleData,
                 targetVariables = c("y", paste0("x", c(1 : 3))),
                 ignoreVariables = "idNum")

