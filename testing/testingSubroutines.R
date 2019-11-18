### Title:    Simulation Subroutines
### Author:   Kyle M. Lang
### Created:  2015-01-01
### Modified: 2019-11-15
### Purpose:  This file contains the subroutines underlying the Monte Carlo
###           simulation to replicate/correct my 2018 ISBA analysis.

###--------------------------------------------------------------------------###
        
## Simulate complete data:
simData <- function(parms, control) {
    nAux <- control$nAux    # Number of auxiliary variables
    nX   <- length(parms$X) # Number of focal predictors

    ## Generate predictor's covariance matrix:
    sigmaZ       <- matrix(parms$zCor, nAux, nAux)
    diag(sigmaZ) <- 1.0
    
    ## Generate auxiliary variables and their regression slopes:
    Z    <- rmvnorm(parms$nObs, rep(0, nAux), sigmaZ)
    zeta <- matrix(parms$zeta, nAux, nX)

    ## Impose sparsity:
    if(control$sparse) zeta[((nAux / 2) + 1) : nAux, ] <- 0.0

    ## Generate focal predictors as X = f(Z):
    s0 <- crossprod(zeta)
    sX <- s0 + diag(((1 / parms$r2$X) - 2) * diag(s0))
    X  <- Z %*% zeta + rmvnorm(parms$nObs, parms$muX, sX)
    
    ## Generate focal outcome as y = f(X, Z):
    beta <- matrix(c(parms$beta, zeta[ , 1]))
    eta  <- cbind(1, X, Z) %*% beta
    sY   <- var(eta) / parms$r2$y - var(eta)
    y    <- eta + rnorm(parms$nObs, 0, sqrt(sY))

    ## Combine data blocks:
    data           <- data.frame(y, X, Z)
    colnames(data) <- c(parms$y, parms$X, paste0("z", 1 : nAux))
    
    data
}

###--------------------------------------------------------------------------###
                                        #data <- missData$data

## Estimate the imputation models:
fitMiModels <- function(data, parms, control, doMice = TRUE) {

    ## Define imputation targets:
    targets <- with(parms, c(y, X))

    if(doMice) {
        ## Run vanilla MICE:
        if(control$expNum < 3)
            naiveMiceOut <- try(
                mice(data      = data,
                     m         = parms$nImps,
                     method    = parms$miceMeth,
                     maxit     = control$miceIters,
                     printFlag = parms$verbose,
                     ridge     = control$miceRidge)
            )

        ## Run MICE using quickpred for variable selection:
        qpPreds <- quickpred(data    = data,
                             mincor  = parms$minPredCor,
                             include = targets)
        
        qpMiceOut <- try(
            mice(data            = data,
                 m               = parms$nImps,
                 method          = parms$miceMeth,
                 predictorMatrix = qpPreds,
                 maxit           = control$miceIters,
                 printFlag       = parms$verbose)
        )

        ## Run MICE using the true imputation model:
        truePreds <- matrix(0,
                            ncol(data),
                            ncol(data),
                            dimnames = list(colnames(data), colnames(data))
                            )

        truePreds[targets, c(targets, control$marPreds)] <- 1
        diag(truePreds)                                  <- 0
        
        trueMiceOut <- try(
            mice(data            = data,
                 m               = parms$nImps,
                 method          = parms$miceMeth,
                 predictorMatrix = truePreds,
                 maxit           = control$miceIters,
                 printFlag       = parms$verbose)
        )

        if(control$expNum < 3) {
            ## Advance the RNG stream:
            .lec.ResetNextSubstream(control$rngStream)
            
            ## Run MIBRR's vanilla MI:
            vanOut <- try(
                vanilla(data        = data,
                        targetVars  = targets,
                        sampleSizes =
                            as.numeric(control$iters[c("finalBurn", "finalGibbs")]),
                        seed        = control$rngStream,
                        userRng     = control$rngStream,
                        verbose     = parms$verbose,
                        ridge       = control$vanRidge)
            )
        }
    }# CLOSE if(doMice)
    
    ## Advance the RNG stream:
    .lec.ResetNextSubstream(control$rngStream)
    
    ## Run MIBEN:
    if(parms$mcem)
        mibenOut <- try(
            miben(data        = data,
                  targetVars  = targets,
                  iterations  = as.numeric(
                      control$iters[c("nMibenEmApprox", "nEmTune")]
                  ),
                  sampleSizes = list(
                      as.numeric(control$iters[c("approxBurn", "approxGibbs")]),
                      as.numeric(control$iters[c("tuneBurn", "tuneGibbs")]),
                      as.numeric(control$iters[c("finalBurn", "finalGibbs")])
                  ),
                  seed        = control$rngStream,
                  userRng     = control$rngStream,
                  verbose     = parms$verbose,
                  control     = list(optCheckKkt    = parms$checkKkt,
                                     optMethod      = parms$optMeth,
                                     optBoundLambda = parms$optBound,
                                     optStrict      = parms$optStrict,
                                     lambda1Starts  = control$lamStarts$miben[[1]],
                                     lambda2Starts  = control$lamStarts$miben[[2]])
                  )
        )
    else
        mibenOut <- try(
            miben(data         = data,
                  targetVars   = targets,
                  sampleSizes  = as.numeric(
                      control$iters[c("finalBurn", "finalGibbs")]
                  ),
                  doMcem       = FALSE,
                  lam1PriorPar = control$l1Par,
                  lam2PriorPar = control$l2Par,
                  seed         = control$rngStream,
                  userRng      = control$rngStream,
                  verbose      = parms$verbose)
        )
    
    ## Advance the RNG stream:
    .lec.ResetNextSubstream(control$rngStream)

    ## Run MIBL:
    if(parms$mcem)
        miblOut <- try(
            mibl(data        = data,
                 targetVars  = targets,
                 iterations  = as.numeric(
                     control$iters[c("nMiblEmApprox", "nEmTune")]
                 ),
                 sampleSizes = list(
                     as.numeric(control$iters[c("approxBurn", "approxGibbs")]),
                     as.numeric(control$iters[c("tuneBurn", "tuneGibbs")]),
                     as.numeric(control$iters[c("finalBurn", "finalGibbs")])
                 ),
                 seed        = control$rngStream,
                 userRng     = control$rngStream,
                 verbose     = parms$verbose,
                 control     = list(lambda1Starts = control$lamStarts$mibl)
                 )
        )
    else
        miblOut <- try(
            mibl(data         = data,
                 targetVars   = targets,
                 sampleSizes  = as.numeric(
                     control$iters[c("finalBurn", "finalGibbs")]
                 ),
                 doMcem       = FALSE,
                 lam1PriorPar = control$l1Par,
                 seed         = control$rngStream,
                 userRng      = control$rngStream,
                 verbose      = parms$verbose)
        )
    
    ## Build and return a combined output object:
    out <- list(miben = mibenOut, mibl = miblOut)
    
    if(doMice) {
        out$quickMice <- qpMiceOut
        out$trueMice  <- trueMiceOut
        
        if(control$expNum < 3) {
            out$naiveMice <- naiveMiceOut
            out$vanilla   <- vanOut
        }
    }
    out
}

###--------------------------------------------------------------------------###

## Extract the multiply imputed datasets:
getImputedData <- function(miOut, parms) {
    if(class(miOut) == "mids") {
        imps <- list()
        for(i in 1 : parms$nImps) imps[[i]] <- complete(miOut, i)
    }
    else {
        imps <- getImpData(miOut, parms$nImps)
    }
    imps
}

###--------------------------------------------------------------------------###

## Fit the analysis models:
fitAnalysisModels <- function(data, y, X) {
    ## Define model formula:
    f <- paste(y, paste0(X, collapse = " + "), sep = " ~ ")

    ## Fit models:
    if(is.data.frame(data)) # Complete data
        out <- lm(f, data = data)
    else {                  # MI data
        tmp <- lapply(data, FUN = function(x) lm(f, data = as.data.frame(x)))
        out <- MIcombine(tmp)
    }
    out
}

###--------------------------------------------------------------------------###

## 'Smart' colMeans function that won't break with vectors:
colMean2 <- function(x, na.rm = FALSE) {
    if(is.vector(x)) mean(x, na.rm = na.rm)
    else             colMeans(x, na.rm = na.rm)
}

###--------------------------------------------------------------------------###

## 'Smart' rowMeans function that won't break with vectors:
rowMean2 <- function(x, na.rm = FALSE) {
    if(is.vector(x)) mean(x, na.rm = na.rm)
    else             rowMeans(x, na.rm = na.rm)
}

###--------------------------------------------------------------------------###

## Compute sufficient statistics:
getStats <- function(data, vars) {
    if(is.data.frame(data)) { # Complete data
        mean <- colMeans(data[ , vars])
        cov  <- cov(data[ , vars])
    }
    else {                    # MI data
        mean <- colMeans(do.call(rbind, data)[ , vars])
        tmp  <- do.call(rbind,
                        lapply(data, function(x, v) cov(x[ , v]), v = vars)
                        )
        tmp  <- aggregate(tmp, by = list(rownames(tmp)), FUN = mean)

        rownames(tmp) <- tmp[ , 1]
        cov           <- tmp[vars, -1]
    }
    rownames(cov) <- vars
    list(mean = mean, cov = cov)
}

###--------------------------------------------------------------------------###

## Check convergence:
converged <- function(obj) class(obj) != "try-error"

###--------------------------------------------------------------------------###

## Get extra info for MIBRR-based models:
getExtras <- function(mibrrOut, parms) {
    out <- list()
    if(converged(mibrrOut)) {
        out$rHats   <- mibrrOut$rHats
        if(parms$mcem) out$lambdas <- mibrrOut$lambdaHistory
        
        check <- parms$mcem && parms$checkKkt & mibrrOut$penalty == 2
        if(check) out$kkt <- mibrrOut$lambdaConv
        
        if(parms$saveParams) {
            out$params <- list()
            for(v in mibrrOut$targetVars)
                out$params[[v]] <- getParams(mibrrOut, v)
        }
    }
    out
}

###--------------------------------------------------------------------------###

## Do the MI analysis phase:
doMiAnalysis <- function(data, parms) {
    if(converged(data)) {# MI converged
        imps  <- getImputedData(miOut = data, parms = parms)
        fit   <- try(fitAnalysisModels(data = imps, y = parms$y, X = parms$X))
        conv  <- c(mi = TRUE, lm = converged(fit))
        stats <- getStats(data = imps, vars = with(parms, c(y, X)))
    }
    else {               # MI crashed
        fit   <- list()
        conv  <- c(mi = FALSE, lm = NA)
        stats <- list(mean = NA, cov = NA)
    }
    list(fit = fit, mean = stats$mean, cov = stats$cov, conv = conv)
}

###--------------------------------------------------------------------------###
                                        #pm <- 0.1

## Execute a single condition of the simulation:
runCondition <- function(pm, compData, control, parms, ...) {
    ## Extract optional parameters:
    args   <- list(...)
    miOnly <- !is.null(args$miOnly) && args$miOnly
    doMice <- !is.null(args$doMice) && args$doMice
    
    ## Incomplete auxiliaries?
    if(parms$pmAux > 0.0) mcarTargets <- paste0("z", 1 : control$nAux)
    else                  mcarTargets <- NA

    ## Impose missingness on the complete data:
    missData <- imposeMissData(data    = compData,
                               targets = list(mar  = with(parms, c(y, X)),
                                              mcar = mcarTargets,
                                              mnar = NA),
                               preds   = control$marPreds,
                               pm      = list(mar  = pm,
                                              mcar = parms$pmAux),
                               snr     = list(mar = parms$marSnr),
                               pattern = parms$marPattern)

                                        #cor(missData$data, use = "pairwise")
                                        #md.pairs(missData$data)$rr / nrow(missData$data)
    
    ## Impute missing data:
    miList <- fitMiModels(data    = missData$data,
                          parms   = parms,
                          control = control,
                          doMice  = doMice)
    
    ## Return early when doing only MI:
    if(miOnly) return(miList)
    
    ## Run MI analysis models:
    resList <- lapply(miList, doMiAnalysis, parms = parms)
    
    ## Extract extra bits:
    extras <- lapply(miList[c("miben", "mibl")], getExtras, parms = parms)
    
    ## Augment the results object:
    resList$miben <- c(resList$miben, extras$miben)
    resList$mibl  <- c(resList$mibl, extras$mibl)

    if(control$expNum < 3) {
        extras          <- getExtras(miList$vanilla, parms)
        resList$vanilla <- c(resList$vanilla, extras)
    }
    
    ## Run complete-data comparison model:
    tmp <- try(fitAnalysisModels(data = compData, y = parms$y, X = parms$X))
    tmp <- list(fit = tmp, conv = converged(tmp))
    
    resList$comp <- c(tmp,
                      getStats(data = compData, vars = with(parms, c(y, X)))
                      )
    
    ## Save the results:
    saveRDS(resList,
            file =
                paste0(parms$outDir,
                       "simRes_", ifelse(control$sparse, "sparse", "dense"),
                       "_v",      control$nAux,
                       "_pm",     100 * pm,
                       "_rep",    control$rp,
                       ".rds")
            )
}# END runCondition()

###--------------------------------------------------------------------------###
rp <- 1

goBabyGo <- function(rp, control, parms) {
    ## Store the current rep index:
    control$rp <- rp

    ## Initialize the L'Ecuyer random number streams:
    .lec.SetPackageSeed(rep(parms$mySeed, 6))
    .lec.CreateStream(paste0(parms$streamStem, c(1 : parms$nReps)))

    ## Associate this replication with it's own stream:
    control$rngStream <- paste0(parms$streamStem, rp)
    .lec.CurrentStream(control$rngStream)

    ## Simulate the complete data:
    compData <- simData(parms = parms, control = control)

    cor(compData)
    colMeans(compData)
    
    ## Apply over percents missing:
    sapply(X        = parms$pmVec,
           FUN      = runCondition,
           compData = compData,
           control  = control,
           parms    = parms,
           doMice   = TRUE)

    ## Clean up the RNG state:
    .lec.CurrentStreamEnd()
    .lec.DeleteStream(paste0(parms$streamStem, 1 : parms$nReps))

    rp # Return the replication number
}# END goBabyGo()

###--------------------------------------------------------------------------###

## Only do MI for a simple level of PM. Used for iteration planning.
iterPlan <- function(rp, pm, parms, control, doMice = TRUE, nChains = 2) {
    ## Store the current rep index:
    control$rp <- rp
    
    ## Initialize the L'Ecuyer random number streams:
    .lec.SetPackageSeed(rep(parms$mySeed, 6))
    .lec.CreateStream(paste0(parms$streamStem, c(1 : parms$nReps)))
    
    ## Associate this replication with it's own stream:
    control$rngStream <- paste0(parms$streamStem, rp)
    .lec.CurrentStream(control$rngStream)
    
    ## Simulate the complete data:
    compData <- simData(parms = parms, control = control)
    
    out <- list()
    for(i in 1 : nChains) {
        ## Perturb starting values:
        control$lamStarts <- lapply(control$lamStarts0,
                                    function(x) x + rnorm(length(x), 0, 0.5 * x)
                                    )
        
        ## Apply over percents missing:
        out[[i]] <- runCondition(pm       = pm,
                                 compData = compData,
                                 control  = control,
                                 parms    = parms,
                                 miOnly   = TRUE,
                                 doMice   = doMice)
    }
    
    ## Clean up the RNG state:
    .lec.CurrentStreamEnd()
    .lec.DeleteStream(paste0(parms$streamStem, 1 : parms$nReps))
    
    out
}

###--------------------------------------------------------------------------###

getRunTime <- function(repTime, nReps, nCores) {
minutes <- ((nReps * repTime) / 60) / nCores
hours   <- minutes / 60
days    <- hours / 24
years   <- days / 365

outMat <- matrix(c(nReps, nCores, minutes, hours, days, years), nrow = 1)
colnames(outMat) <- c("Reps", "Cores", "Minutes", "Hours", "Days", "Years")
outMat
}

###--------------------------------------------------------------------------###

## Broadcast the library function of a list of packages:
applyLib <- function(pkgList)
    lapply(pkgList, library, character.only = TRUE, logical = TRUE)

###--------------------------------------------------------------------------###

justVanilla <- function(rp, pm, parms, control) {
    ## Store the current rep index:
    control$rp <- rp
    
    ## Initialize the L'Ecuyer random number streams:
    .lec.SetPackageSeed(rep(parms$mySeed, 6))
    .lec.CreateStream(paste0(parms$streamStem, c(1 : parms$nReps)))
    
    ## Associate this replication with it's own stream:
    control$rngStream <- paste0(parms$streamStem, rp)
    .lec.CurrentStream(control$rngStream)

    ## Simulate some complete data:
    compData <- simData(parms, control)
    
    ## Impose missingness on the complete data:
    missData <- imposeMissData(data    = compData,
                               targets = list(mar  = with(parms, c(y, X)),
                                              mcar = NA,
                                              mnar = NA),
                               preds   = control$marPreds,
                               pm      = list(mar  = pm,
                                              mcar = parms$pmAux),
                               snr     = list(mar = parms$marSnr),
                               pattern = parms$marPattern)$data
    
    ## Run vanilla MI:
    vanOut <- try(
        vanilla(data        = missData,
                targetVars  = with(parms, c(y, X)),
                sampleSizes =
                    as.numeric(control$iters[c("finalBurn", "finalGibbs")]),
                seed        = control$rngStream,
                userRng     = control$rngStream,
                verbose     = parms$verbose,
                ridge       = control$vanRidge)
    )
}

###--------------------------------------------------------------------------###

plotTrace <- function(x, y, what) {
    x <- x[[tolower(what)]]
    y <- y[[tolower(what)]]
    
    if(tolower(what) == "sigma") {
        yLim <- range(x, y)
        plot(x    = x,
             ylim = yLim,
             type = "l",
             col  = "red",
             ylab = i,
             main = paste0("Trace of ", what)
             )
        lines(y, col = "blue")
    }
    else {
        for(v in 1 : ncol(x)) {
            yLim <- range(x[ , v], y[ , v])
            plot(x    = x[ , v],
                 ylim = yLim,
                 type = "l",
                 col  = "red",
                 ylab = what,
                 main = switch(tolower(what),
                               lambda = paste0("Trace of ", what, v),
                               paste0("Trace of ",
                                      what,
                                      " for ",
                                      colnames(x)[v])
                               )
                 )
            lines(y[ , v], col = "blue")
            
            readline("Press any key to continue. ")
        }
    }
}

###--------------------------------------------------------------------------###

                                        #library(SURF)
                                        #?simCovData
                                        #?imposeMissData

## Simulate some very simple data:
genSimpleData <- function(n, cov, pm) {
    ## Generate complete data:
    dat0 <- SURF::simCovData(nObs = n, sigma = cov, nVars = 3)

    ## Impose missing:
    r1 <- dat0$x3 < quantile(dat0$x3, pm)
    r2 <- dat0$x3 > quantile(dat0$x3, (1 - pm))
    
    dat0[r1, "x1"] <- NA
    dat0[r2, "x2"] <- NA

    dat0
}

###--------------------------------------------------------------------------###

## Estimate the imputation models:
runMibrr <- function(data, parms) {

    ## Define imputation targets:
    targets <- colnames(data)[colMeans(is.na(data)) > 0]
    
    ## Advance the RNG stream:
    .lec.ResetNextSubstream(parms$rngStream)
    
    ## Run MIBEN:
    mibenOut <- try(
        miben(data        = data,
              targetVars  = targets,
              iterations  = as.numeric(
                  parms$iters[c("nMibenEmApprox", "nEmTune")]
              ),
              sampleSizes = list(
                  as.numeric(parms$iters[c("approxBurn", "approxGibbs")]),
                  as.numeric(parms$iters[c("tuneBurn", "tuneGibbs")]),
                  as.numeric(parms$iters[c("finalBurn", "finalGibbs")])
              ),
              seed        = parms$rngStream,
              userRng     = parms$rngStream,
              verbose     = parms$verbose,
              parms       = list(optCheckKkt    = parms$checkKkt,
                                 optMethod      = parms$optMeth,
                                 optBoundLambda = parms$optBound,
                                 optStrict      = parms$optStrict,
                                 lambda1Starts  = parms$lamStarts$miben[[1]],
                                 lambda2Starts  = parms$lamStarts$miben[[2]])
              )
    )
    
    ## Advance the RNG stream:
    .lec.ResetNextSubstream(parms$rngStream)

    ## Run MIBL:
    miblOut <- try(
        mibl(data        = data,
             targetVars  = targets,
             iterations  = as.numeric(
                 parms$iters[c("nMiblEmApprox", "nEmTune")]
             ),
             sampleSizes = list(
                 as.numeric(parms$iters[c("approxBurn", "approxGibbs")]),
                 as.numeric(parms$iters[c("tuneBurn", "tuneGibbs")]),
                 as.numeric(parms$iters[c("finalBurn", "finalGibbs")])
             ),
             seed        = parms$rngStream,
             userRng     = parms$rngStream,
             verbose     = parms$verbose,
             parms       = list(lambda1Starts = parms$lamStarts$mibl)
             )
    )
    
    ## Return a combined output object:
    list(miben = mibenOut, mibl = miblOut)
}

###--------------------------------------------------------------------------###

## Only do MI for a simple level of PM. Used for iteration planning.
testMcem <- function(rp, pm, parms, nChains = 2) {
    ## Store the current rep index:
    parms$rp <- rp
    
    ## Initialize the L'Ecuyer random number streams:
    .lec.SetPackageSeed(rep(parms$mySeed, 6))
    .lec.CreateStream(paste0(parms$streamStem, c(1 : parms$nReps)))
    
    ## Associate this replication with it's own stream:
    parms$rngStream <- paste0(parms$streamStem, rp)
    .lec.CurrentStream(parms$rngStream)
    
    ## Simulate the complete data:
    dat0 <- genSimpleData(n = parms$nObs, cov = parms$cx, pm = pm)
    
    out <- list()
    for(i in 1 : nChains) {
        ## Perturb starting values:
        parms$lamStarts <- lapply(parms$lamStarts0,
                                  function(x) x + rnorm(length(x), 0, 0.5 * x)
                                  )
        
        ## Apply over percents missing:
        out[[i]] <- runMibrr(data = dat0, parms = parms)
    }
    
    ## Clean up the RNG state:
    .lec.CurrentStreamEnd()
    .lec.DeleteStream(paste0(parms$streamStem, 1 : parms$nReps))
    
    out
}

###--------------------------------------------------------------------------###
