### Title:    Build R Packages
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2018-MAY-02
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))

library(RcppEigen)

ver <- "0.0.0.9005-0000"

system("rm source/MIBRR/src/RcppExports.cpp \
        rm source/MIBRR/R/RcppExports.R")
Rcpp::compileAttributes("source/MIBRR")

system("R CMD build source/MIBRR")
install.packages(paste0("MIBRR_", ver, ".tar.gz"),
                 repos = NULL,
                 type  = "source")

system(paste0("mv MIBRR_", ver, ".tar.gz builds/"))
