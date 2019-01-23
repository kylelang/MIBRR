### Title:    Build R Packages
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2019-JAN-23
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))

library(RcppEigen)

ver <- "0.3.0.9000"

## Clean up:
system("rm source/MIBRR/src/RcppExports.cpp \
        rm source/MIBRR/R/RcppExports.R")

## Generate the RccpExports.cpp file:
Rcpp::compileAttributes("source/MIBRR")

## Build the MIBRR package:
system("R CMD build source/MIBRR")

## Run CRAN checks:
system(paste0("R CMD check MIBRR_", ver, ".tar.gz"))

## Install the MIBRR package:
install.packages(paste0("MIBRR_", ver, ".tar.gz"),
                 repos = NULL,
                 type  = "source")

## Move the tar-ball to the builds directory:
system(paste0("mv MIBRR_", ver, ".tar.gz builds/"))

## Generate a PDF version of the manual:
system(paste0("R CMD Rd2pdf --force --output='./documentation/MIBRR.pdf' ",
             .libPaths()[1],
             "/MIBRR")
       )
