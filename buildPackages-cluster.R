### Title:    Build/Install the MIBRR R Package
### Author:   Kyle M. Lang
### Created:  2014-12-07
### Modified: 2019-12-06

rm(list = ls(all = TRUE))

library(RcppEigen)

ver  <- "0.3.3.9001"
prod <- FALSE

## Clean up:
system("rm source/MIBRR/src/RcppExports.cpp \
        rm source/MIBRR/R/RcppExports.R")

## Generate the RccpExports.cpp file:
Rcpp::compileAttributes("source/MIBRR")

## Build the MIBRR package:
system("R CMD build source/MIBRR")

## Move the tar-ball to the builds directory:
system(paste0("scp MIBRR_", ver, ".tar.gz kmlang@lisa.surfsara.nl:/home/kmlang/rPackages"))
