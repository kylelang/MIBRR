### Title:    Build R Packages
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2017-SEP-30
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))

install.packages(c("optimx"), repos = "http://cloud.r-project.org")

library(RcppEigen)

system("rm source/mibrr/src/RcppExports.cpp \
        rm source/mibrr/R/RcppExports.R")
Rcpp::compileAttributes("source/mibrr")
system("R CMD build source/mibrr")
install.packages("mibrr_0.0.0.9001.tar.gz", repos = NULL, type = "source")
