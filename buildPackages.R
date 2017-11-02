### Title:    Build R Packages
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2017-OCT-27
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))

library(RcppEigen)

system("rm source/mibrr/src/RcppExports.cpp \
        rm source/mibrr/R/RcppExports.R")
Rcpp::compileAttributes("source/mibrr")
system("R CMD build source/mibrr")
install.packages("mibrr_0.0.0.9005.tar.gz", repos = NULL, type = "source")
