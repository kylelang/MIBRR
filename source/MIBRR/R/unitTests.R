### Title:    R-Based Unit Tests for MIBRR
### Author:   Kyle M. Lang
### Created:  2010-01-23
### Modified: 2019-12-12

##--------------------- COPYRIGHT & LICENSING INFORMATION --------------------##
##  Copyright (C) 2019 Kyle M. Lang <k.m.lang@uvt.nl>                         ##
##                                                                            ##
##  This file is part of MIBRR.                                               ##
##                                                                            ##
##  This program is free software: you can redistribute it and/or modify it   ##
##  under the terms of the GNU General Public License as published by the     ##
##  Free Software Foundation, either version 3 of the License, or (at you     ##
##  option) any later version.                                                ##
##                                                                            ##
##  This program is distributed in the hope that it will be useful, but       ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of                ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General  ##
##  Public License for more details.                                          ##
##                                                                            ##
##  You should have received a copy of the GNU General Public License along   ##
##  with this program. If not, see <http://www.gnu.org/licenses/>.            ##
##----------------------------------------------------------------------------##

### Random Number Samplers ###

## Compare MIBRR's multivariate normal sampler to a reference implementation
## using the Kolmogorov-Smirnov Statistic:
testMvn <- function(n = 1000, seed = NA, ...) {
    pars <- list(...)[[1]]
      
    if(!is.null(pars$mu)) mu <- pars$mu
    else                  mu <- rep(0, 2)
    
    if(!is.null(pars$sigma)) sigma <- pars$sigma
    else                     sigma <- diag(2)
    
    if(is.na(seed)) seed <- floor(runif(1, 1e5, 1e6))
    
    sam1 <- drawMvn(n, mu, sigma, seed)
    sam2 <- rmvnorm(n, mu, sigma)
    
    out <- list()
    for(v in 1 : length(mu))
        out[[v]] <- ks.test(x = sam1[ , v], y = sam2[ , v])
    
    out
}

###--------------------------------------------------------------------------###

## Compare MIBRR's inverse gamma sampler to a reference implementation using the
## Kolmogorov-Smirnov Statistic:
testInvGamma <- function(n = 1000, seed = NA, ...) {
    pars <- list(...)[[1]]

    if(!is.null(pars$shape)) shape <- pars$shape
    else                     shape <- 1
    
    if(!is.null(pars$scale)) scale <- pars$scale
    else                     scale <- 1

    if(is.na(seed)) seed <- floor(runif(1, 1e5, 1e6))
    
    sam1 <- drawInvGamma(n, shape, scale, seed)
    sam2 <- rinvgamma(n, shape, scale)
    
    ks.test(x = sam1, y = sam2)
}

###--------------------------------------------------------------------------###

## Compare MIBRR's scaled inverse chi-squared sampler to a reference
## implementation using the Kolmogorov-Smirnov Statistic:
testInvChiSq <- function(n = 1000, seed = NA, ...) {
    pars <- list(...)[[1]]
    
    if(!is.null(pars$df)) df <- pars$df
    else                  df <- 1
    
    if(!is.null(pars$scale)) scale <- pars$scale
    else                     scale <- 1

    if(is.na(seed)) seed <- floor(runif(1, 1e5, 1e6))
    
    sam1 <- drawInvChiSq(n, df, scale, seed)
    sam2 <- rinvchisq(n, df, scale)
    
    ks.test(x = sam1, y = sam2)
}

###--------------------------------------------------------------------------###

## Compare MIBRR's inverse Gaussian sampler to a reference implementation using
## the Kolmogorov-Smirnov Statistic:
testInvGauss <- function(n = 1000, seed = NA, ...) {
    pars <- list(...)[[1]]
    
    if(!is.null(pars$mu)) mu <- pars$mu
    else                  mu <- 1
    
    if(!is.null(pars$lambda)) lambda <- pars$lambda
    else                      lambda <- 1
    
    if(is.na(seed)) seed <- floor(runif(1, 1e5, 1e6))
    
    sam1 <- drawInvGauss(n, mu, lambda, seed)
    sam2 <- rinvgaussian(n, mu, lambda)
    
    ks.test(x = sam1, y = sam2)
}

###--------------------------------------------------------------------------###

## Compare MIBRR's generalized inverse Gaussian sampler to a reference
## implementation using the Kolmogorov-Smirnov Statistic:
testGig <- function(n = 1000, seed = NA, ...) {
    pars <- list(...)[[1]]
    
    if(!is.null(pars$lambda)) lambda <- pars$lambda
    else                      lambda <- 1

    if(!is.null(pars$chi)) chi <- pars$chi
    else                   chi <- 1
    
    if(!is.null(pars$psi)) psi <- pars$psi
    else                   psi <- 1

    if(is.na(seed)) seed <- floor(runif(1, 1e5, 1e6))
    
    sam1 <- drawGig(n, lambda, chi, psi, seed)
    sam2 <- rgig(n, c(lambda, chi, psi))

    ks.test(x = sam1, y = sam2)
}

###--------------------------------------------------------------------------###

## Run each sampler test once:
samTestUnit <- function(n = 1000, ...) {
    pars <- list(...)
    
    ## Generate an RNG seed:
    seed <- floor(runif(1, 1e5, 1e6))
    
    mvn        <- sapply(testMvn(n, seed, pars$mvn), "[[", x = "statistic")
    names(mvn) <- paste0("mvn", 1 : length(mvn))

    out <- c(mvn,
             invGamma = testInvGamma(n, seed, pars$invGamma)$statistic,
             invGauss = testInvGauss(n, seed, pars$invGauss)$statistic,
             gig      = testGig(n, seed, pars$gig)$statistic,
             invChiSq = testInvChiSq(n, seed, pars$invChiSq)$statistic)

    names(out) <- gsub("\\.D", "", names(out))
    out
}

###--------------------------------------------------------------------------###

## Run a Monte Carlo test of the samplers:
testSamplers <- function(reps = 1000, n = 1000, ...) {
    out <- list()
    for(i in 1 : reps) out[[i]] <- samTestUnit(n, ...)

    note <- strwrap("The following table summarizes the Monte Carlo sampling distributions of KS statistics that compare samples generated by MIBRR's samplers to those generated by R-based reference implementations from the mvtnorm, LaplacesDemon, and HyperbolicDist packages.",
                    width = 81)
    
    for(i in note) message(i)
    
    message(paste0("\nNo. of Monte Carlo replications = ", reps, "\n"))
    message(paste0("No. of observations in each sample = ", n, "\n"))
    
    out <- do.call(rbind, out)
    out <- rbind(colMeans(out),
                 apply(out, 2, sd),
                 apply(out, 2, quantile, probs = c(0.9, 0.95, 0.99, 0.999)),
                 apply(out, 2, max)
                 )
    
    rownames(out) <- c("Mean",
                       "SD",
                       paste0(c(90, 95, 99, 99.9), "th %ile"),
                       "Maximum")
    out
}

###--------------------------------------------------------------------------###

### Missing Data Indexing ###

testMissIndex <- function(data) {
    ## Summarize missing data in the raw data file:
    missList0   <- lapply(data, function(x) which(is.na(x)) - 1)
    respCounts0 <- colSums(!is.na(data))
    
    ## Initialize a MibrrFit object:
    tmp <- suppressWarnings(miben(data = data, initOnly = TRUE))
    
    respCounts <- nrow(tmp$data) - tmp$missCounts

    ## Check R-level missing data accounting:
    check <- all.equal(tmp$missList, missList0[names(tmp$missList)])
    if(!check) stop("R-level missing data indexing is broken.")
    
    check <- all(respCounts == respCounts0[names(respCounts)])
    if(!check) stop("R-level response counting is broken.")
    
    ## Check C++-level missing data accounting:
    for(v in 1 : ncol(tmp$data)) {
        poiOut <- printObsIndices(data        = as.matrix(tmp$data),
                                  missIndices = tmp$missList,
                                  respCounts  = respCounts,
                                  targetIndex = v - 1)
        
        pmiOut <- printMissIndices(data        = as.matrix(tmp$data),
                                   missIndices = tmp$missList,
                                   respCounts  = respCounts,
                                   targetIndex = v - 1)
        
        ## Any overlap?
        t1 <- intersect(pmiOut, poiOut)
        if(length(t1) > 0) stop("Miss and obs indices are not disjoint")
        
        ## All rows maintained?
        outInds <- c(pmiOut, poiOut)
        t2 <- setdiff(c(1 : nrow(tmp$data)), outInds + 1)
        if(length(t2) > 0) stop("Some rows dropped")
        
        ## In == Out?
        t3 <- setdiff(pmiOut, tmp$missList[[v]])
        if(length(t3) > 0) stop("Miss indices broken")
    }# END for(v in 1 : ncol(tmp$data))

    ## Success; no errors
    TRUE
}

###--------------------------------------------------------------------------###

### Data Manipulation ###

testDataProcessing <- function(data) {
    ## Initialize a MibrrFit object:
    obj <- suppressWarnings(miben(data = data, initOnly = TRUE))

    dat1       <- obj$data
    missList   <- obj$missList
    respCounts <- nrow(dat1) - obj$missCounts
    
    ## Check initial data integrity:
    diff <- sum(dat1 - data, na.rm = TRUE)
    if(diff != 0) stop("MibrrFit$data doesn't match raw data.")

    ## Check data subsetting:
    for(v in 1 : ncol(data)) {
        y0  <- data[ , v]
        X0o <- data[!is.na(y0) , -v] 
        X0m <- data[is.na(y0), -v]

        ## Extract "training set" predictors:
        X1o <- getX(data        = as.matrix(dat1),
                    missIndices = missList,
                    respCounts  = respCounts,
                    noMiss      = FALSE,
                    xOnly       = TRUE,
                    obsY        = TRUE,
                    scale       = FALSE,
                    targetIndex = v - 1)
        
        diff <- sum(X1o - X0o, na.rm = TRUE)
        if(diff != 0) stop("X subsetting is broken for observed y.")

        ## Extract "testing set" predictors:
        X1m <- getX(data        = as.matrix(dat1),
                    missIndices = missList,
                    respCounts  = respCounts,
                    noMiss      = FALSE,
                    xOnly       = TRUE,
                    obsY        = FALSE,
                    scale       = FALSE,
                    targetIndex = v - 1)
        
        diff <- sum(X1m - X0m, na.rm = TRUE)
        if(diff != 0) stop("X subsetting is broken for missing y.")
        
        y1 <- getY(data        = as.matrix(dat1),
                   missIndices = missList,
                   respCounts  = respCounts,
                   noMiss      = FALSE,
                   targetIndex = v - 1)

        diff <- sum(y1 - na.omit(y0))
        if(diff != 0) stop("y subsetting is broken.")
    }

    ## Success; no errors:
    TRUE
}

###--------------------------------------------------------------------------###

### Data Scaling ###

testDataScaling <- function(data) {
    ## Initialize a MibrrFit object:
    obj <- suppressWarnings(miben(data = data, initOnly = TRUE))

    dat1       <- obj$data
    missList   <- obj$missList
    respCounts <- nrow(dat1) - obj$missCounts
    
    ## Check data subsetting:
    for(v in 1 : ncol(dat1)) {
        ## Extract "training set" predictors:
        Xo <- getX(data        = as.matrix(dat1),
                   missIndices = missList,
                   respCounts  = respCounts,
                   noMiss      = FALSE,
                   xOnly       = TRUE,
                   obsY        = TRUE,
                   scale       = TRUE,
                   targetIndex = v - 1)

        m1    <- colMeans(Xo)
        check <- all.equal(m1, rep(0, ncol(dat1) - 1))

        if(!check) stop("Training set mean centering is broken.")
        
        s1    <- apply(Xo, 2, sd)
        check <- all.equal(s1, rep(1, ncol(dat1) - 1))

        if(!check) stop("Training set scaling is broken.")
        
        ## Extract scaled "testing set" predictors:
        Xm <- getX(data        = as.matrix(dat1),
                   missIndices = missList,
                   respCounts  = respCounts,
                   noMiss      = FALSE,
                   xOnly       = TRUE,
                   obsY        = FALSE,
                   scale       = TRUE,
                   targetIndex = v - 1)

        ## Compute centers and scales:
        tmp <- getX(data        = as.matrix(dat1),
                    missIndices = missList,
                    respCounts  = respCounts,
                    noMiss      = FALSE,
                    xOnly       = TRUE,
                    obsY        = TRUE,
                    scale       = FALSE,
                    targetIndex = v - 1)
        
        m0 <- matrix(colMeans(tmp), nrow(Xm), ncol(Xm), byrow = TRUE)
        s0 <- matrix(apply(tmp, 2, sd), nrow(Xm), ncol(Xm), byrow = TRUE)

        ## Revert scaling of the testing set predictors:
        Xm0 <- Xm * s0 + m0

        ## Extract raw "testing set" predictors:
        Xm1 <- getX(data        = as.matrix(dat1),
                    missIndices = missList,
                    respCounts  = respCounts,
                    noMiss      = FALSE,
                    xOnly       = TRUE,
                    obsY        = FALSE,
                    scale       = FALSE,
                    targetIndex = v - 1)

        ## Check that we can recover the raw data:
        check <- all.equal(Xm1, Xm0)
        if(!check) stop("Testing set scaling is broken.")
    }
    
    ## Success; no errors:
    TRUE
}

###--------------------------------------------------------------------------###

### Missing Data Filling ###

testMissFill <- function(data) {
    ## Initialize a MibrrFit object:
    obj <- suppressWarnings(miben(data = data, initOnly = TRUE))

    missList   <- obj$missList
    respCounts <- nrow(data) - obj$missCounts

    ## Check C++-level missing data replacement ##

    for(v in 1 : ncol(data)) {
        ## Create some dummy imputations:
        imps <- runif(obj$missCounts[v])

        ## Manually fill missing:
        y0            <- data[ , v]
        y0[is.na(y0)] <- imps

        ## Fill the missing in C++ code:
        y1 <- printFilledY(imps,
                           as.matrix(obj$data),
                           missList,
                           respCounts,
                           v - 1)
        
        ## Compare manual and C++ versions:
        diff <- sum(y1 - y0)
        if(diff != 0) stop("C++-level missing value replacement is broken.")
    }
    
    ## Check R-level missing data replacement ##
    
    ## Generate a full MibrrFit object:
    obj <- suppressWarnings(
        miben(data         = data,
              sampleSizes  = c(5, 5), 
              doMcem       = FALSE,
              lam1PriorPar = c(1.0, 0.1),
              lam2PriorPar = c(1.0, 0.1),
              control      = list(checkConv = FALSE),
              verbose      = FALSE)
    )
    
    for(rep in 1 : 3) {
        dat0 <- data # Raw data
        for(v in colnames(obj$data)) {
            nMiss <- obj$missCounts[v]
            
            ## Generate some dummy imputations:
            imp0 <- runif(nMiss)
            
            ## Replace MibrrFit imps with dummies:
            pars      <- obj$chains[[1]]$parameters[[v]]
            pars$imps <- matrix(imp0,
                                length(pars$sigma),
                                nMiss,
                                byrow = TRUE)
            
            ## Replace raw data's missing with dummy imputations:
            dat0[is.na(dat0[ , v]), v] <- imp0
        }
        
        ## Get an imputed dataset from the MibrrFit object:
        dat1 <- obj$getImpDataset()
        
        ## Compare the manual and MibrrFit imputations:
        check <- all(dat1 == dat0)
        if(!check) stop("R-level missing data replacement is broken.")
    }
    
    ## Success, no errors:
    TRUE
}
