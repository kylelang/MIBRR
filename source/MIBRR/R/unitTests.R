### Title:    R-Based Unit Tests for MIBRR
### Author:   Kyle M. Lang
### Created:  2010-JAN-23
### Modified: 2019-JAN-24

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

## MVN sampler:
testMvn <- function(n = 1000, p = 0.1, seed = NA) {
    mu          <- rep(10, 2)
    sigma       <- matrix(5, 2, 2)
    diag(sigma) <- 20

    if(is.na(seed)) seed <- floor(runif(1, 1e5, 1e6))

    sam1 <- drawMvn(n, mu, sigma, seed)
    sam2 <- rmvnorm(n, mu, sigma)
    
    pVec <- rep(NA, 2)
    for(v in 1 : 2)
        pVec[v] <- ks.test(x = sam1[ , v], y = sam2[ , v])$p.value
    
    all(pVec >= p)
}

###--------------------------------------------------------------------------###

## Inverse gamma sampler:
testInvGamma <- function(n = 1000, p = 0.1, seed = NA) {
    shape <- 10
    scale <- 10
    
    if(is.na(seed)) seed <- floor(runif(1, 1e5, 1e6))

    sam1 <- drawInvGamma(n, shape, scale, seed)
    sam2 <- rinvgamma(n, shape, scale)
    
    ks.test(x = sam1, y = sam2)$p.value >= p
}

###--------------------------------------------------------------------------###

## Inverse Gaussian sampler:
testInvGauss <- function(n = 1000, p = 0.1, seed = NA) {
    mu  <- 1
    lam <- 2

    if(is.na(seed)) seed <- floor(runif(1, 1e5, 1e6))
    
    sam1 <- drawInvGauss(n, mu, lam, seed)
    sam2 <- rinvgaussian(n, mu, lam)
    
    ks.test(x = sam1, y = sam2)$p.value >= p
}

###--------------------------------------------------------------------------###

## Generalized inverse Gaussian sampler:
testGig <- function(n = 1000, p = 0.1, seed = NA) {
    lam <- 1
    chi <- 2
    psi <- 2

    if(is.na(seed)) seed <- floor(runif(1, 1e5, 1e6))
    
    sam1 <- drawGig(n, lam, chi, psi, seed)
    sam2 <- rgig(n, c(lam, chi, psi))

    ks.test(x = sam1, y = sam2)$p.value >= p
}

###--------------------------------------------------------------------------###

## Scaled inverse chi-squared sampler:
testScaledInvChiSq <- function(n = 1000, p = 0.1, seed = NA) {
    df    <- 100
    scale <- 10
 
    if(is.na(seed)) seed <- floor(runif(1, 1e5, 1e6))
    
    sam1 <- drawScaledInvChiSq(n, df, scale, seed)
    sam2 <- rinvchisq(n, df, scale)
    
    ks.test(x = sam1, y = sam2)$p.value >= p
}

###--------------------------------------------------------------------------###

## Incomplete gamma calculation:
testIncGamma <- function() {
    shape <- 10
    cut   <- 5
    
    out1 <- calcIncGamma(shape, cut, FALSE)
    out2 <- pgamma(q = cut, shape = shape, lower = FALSE) * gamma(shape)
    
    all.equal(out1, out2)
}

###--------------------------------------------------------------------------###

## Run a single iteration of all sampler tests:
combSamTests <- function(n = 1000, p = 0.1, seed = NA) {
    c(mvn       = testMvn(n, p, seed),
      invGamma  = testInvGamma(n, p, seed),
      invGauss  = testInvGauss(n, p, seed),
      gig       = testGig(n, p, seed),
      sInvChiSq = testScaledInvChiSq(n, p, seed),
      incGamma  = testIncGamma()
      )
}

###--------------------------------------------------------------------------###

## Run a Monte Carlo test of the samplers:
testSamplers <- function(reps = 5000, n = 1000, p = 0.1, seed = NA) {
    out <- list()
    for(i in 1 : reps) out[[i]] <- combSamTests(n, p, seed)
    
    res <- rbind(colMeans(do.call(rbind, out)),
                 c(1 - 2*p, rep(1 - p, 4), 1)
                 )
    
    rownames(res) <- c("Estimated", "Nominal")
    res
}

###--------------------------------------------------------------------------###

### Missing Data Indexing ###

testMissIndex <- function() {
    data(mibrrExampleData)
    
    ## Summarize missing data in the raw data file:
    missList0   <- lapply(mibrrExampleData, function(x) which(is.na(x)) - 1)
    respCounts0 <- colSums(!is.na(mibrrExampleData))
    
    ## Initialize a MibrrFit object:
    tmp <- miben(data       = mibrrExampleData,
                 iterations = c(30, 10),
                 targetVars = c("y", paste0("x", c(1 : 3))),
                 ignoreVars = "idNum",
                 initOnly   = TRUE)
    
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

testDataProcessing <- function() {
    data(mibrrExampleData)
    dat0 <- mibrrExampleData[ , -1]
    
    ## Initialize a MibrrFit object:
    obj <- miben(data       = mibrrExampleData,
                 targetVars = c("y", paste0("x", c(1 : 3))),
                 ignoreVars = "idNum",
                 initOnly   = TRUE)

    dat1       <- obj$data
    missList   <- obj$missList
    respCounts <- nrow(dat1) - obj$missCounts

    ## Check initial data integrity:
    diff <- sum(dat1 - dat0, na.rm = TRUE)
    if(diff != 0) stop("MibrrFit$data doesn't match raw data.")

    ## Check data subsetting:
    for(v in 1 : ncol(dat0)) {
        y0  <- dat0[ , v]
        X0o <- dat0[!is.na(y0) , -v] 
        X0m <- dat0[is.na(y0), -v]

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

testDataScaling <- function() {
    data(mibrrExampleData)
       
    ## Initialize a MibrrFit object:
    obj <- miben(data       = mibrrExampleData,
                 targetVars = c("y", paste0("x", c(1 : 3))),
                 ignoreVars = "idNum",
                 initOnly   = TRUE)

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
