### Title:    MibrrFit Reference Class Definition
### Author:   Kyle M. Lang
### Created:  2017-NOV-28
### Modified: 2019-JAN-31
### Note:     MibrrFit is the metadata class for the MIBRR package

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

MibrrFit <- setRefClass("MibrrFit",
                        fields = list(
                            data              = "data.frame",
                            targetVars        = "character",
                            ignoreVars        = "character",
                            iterations        = "integer",
                            sampleSizes       = "list",
                            missCode          = "integer",
                            seed              = "ANY",
                            doImp             = "logical",
                            doMcem            = "logical",
                            checkConv         = "logical",
                            verbose           = "logical",
                            convThresh        = "numeric",
                            lambda1Starts     = "numeric",
                            lambda2Starts     = "numeric",
                            l1Pars            = "numeric",
                            l2Pars            = "numeric",
                            usePcStarts       = "logical",
                            smoothingWindow   = "integer",
                                        #scale             = "logical",
                            minPredCor        = "numeric",
                            miceIters         = "integer",
                            miceRidge         = "numeric",
                            miceMethod        = "character",
                            preserveStructure = "logical",
                            optTraceLevel     = "integer",
                            optCheckKkt       = "logical",
                            optMethod         = "character",
                            optBoundLambda    = "logical",
                            gibbsOut          = "list",
                            ignoredColumns    = "data.frame",
                            rawNames          = "character",
                            impRowsPool       = "integer",
                            missList          = "list",
                            nChains           = "integer",
                            rHats             = "list",
                            lambdaMat         = "matrix",
                            lambdaHistory     = "list",
                            betaStarts        = "matrix",
                            tauStarts         = "matrix",
                            sigmaStarts       = "numeric",
                            userMissCode      = "logical",
                            missCounts        = "integer",
                            nTargets          = "integer",
                            nVar              = "integer",
                            nPreds            = "integer",
                            nObs              = "integer",
                            totalIters        = "integer",
                            rng0              = "character",
                            userRng           = "character",
                            ridge             = "numeric",
                            penalty           = "integer",
                            savePpSams        = "logical"
                        )
                        )


MibrrFit$methods(
             
################################ CONSTRUCTOR ###################################
             
             initialize =
                 function(data        = data.frame(NULL),
                          targetVars  = "",
                          ignoreVars  = "",
                          iterations  = as.integer(NA),
                          sampleSizes = list(),
                          missCode    = as.integer(NA),
                          doImp       = as.logical(NA),
                          doMcem      = as.logical(NA),
                          verbose     = as.logical(NA),
                          seed        = NULL,
                          userRng     = "",
                          ridge       = 0.0,
                          penalty     = 2L)
                 {
                     "Initialize an object of class MibrrFit"
                     data              <<- data
                     targetVars        <<- targetVars
                     ignoreVars        <<- ignoreVars
                     iterations        <<- iterations
                     sampleSizes       <<- sampleSizes
                     missCode          <<- missCode
                     doImp             <<- doImp
                     doMcem            <<- doMcem
                     checkConv         <<- TRUE
                     verbose           <<- verbose
                     convThresh        <<- 1.1
                     usePcStarts       <<- FALSE
                                        #scale             <<- TRUE
                     minPredCor        <<- 0.3
                     miceIters         <<- 10L
                     miceRidge         <<- 1e-4
                     miceMethod        <<- "pmm"
                     preserveStructure <<- TRUE
                     optTraceLevel     <<- 0L
                     optCheckKkt       <<- TRUE
                     optMethod         <<- "L-BFGS-B"
                     optBoundLambda    <<- TRUE
                     nChains           <<- 1L
                     seed              <<- seed
                     userRng           <<- userRng
                     ridge             <<- ridge
                     penalty           <<- penalty
                     savePpSams        <<- FALSE
                 },
             
################################### MUTATORS ###################################
             
             setData = function(dat1, vars = NULL) {
                 if(is.null(vars))
                     data[ , colnames(dat1)] <<- dat1
                 else
                     data[ , vars] <<- dat1[ , vars]
             },
             
###--------------------------------------------------------------------------###
             
             setControl = function(x) {
                 "Assign the control parameters"
                 ints <- c("smoothingWindow", "miceIters", "optTraceLevel")
                 
                 for(n in names(x)) {
                     if(n %in% ints) field(n, as.integer(x[[n]]))
                     else            field(n, x[[n]])
                 }
             },

###--------------------------------------------------------------------------###
             
             setLambdaParams = function(l1, l2) {
                 l1Pars <<- l1
                 l2Pars <<- l2
             },
             
###--------------------------------------------------------------------------###

             saveSeed = function(seed0) {
                 if(is.character(seed0))
                     seed <<- list(name = seed0, value = .lec.GetState(seed0))
                 else
                     seed <<- list(name = NULL, value = seed0)
             },

###--------------------------------------------------------------------------###

             setupRng = function() {
                 "Start up the l'ecuyer RNG streams"
                 if(is.null(seed))           # No user-supplied seed
                     saveSeed(round(runif(6) * 1e8))
                 else if(is.numeric(seed))   # Numeric user-supplied seed
                     saveSeed(rep(seed, length.out = 6))
                 else if(is.character(seed)) # 'seed' names an extant RNG stream
                     saveSeed(seed)
                 else
                     stop(paste0("The value provided for 'seed' has the wrong type (i.e., ",
                                 class(seed),
                                 ").")
                          )
                 
                 ## Stop generating random numbers from the user's stream:
                 check <- length(userRng) > 0 & userRng != ""
                 if(check) .lec.CurrentStreamEnd()
                 
                 ## Seed the l'ecuyer RNG:
                 .lec.SetPackageSeed(seed$value)
                 
                 ## Generate the l'ecuyer RNG streams:
                 .lec.CreateStream(
                     paste0("mibrrStream", c(0 : length(targetVars)))
                 )

                 ## Set the 'master' RNG stream:
                 rng0 <<- .lec.CurrentStream("mibrrStream0")
             },

###--------------------------------------------------------------------------###

             cleanRng = function() {
                 "Clean up the RNG state"
                 .lec.CurrentStreamEnd(rng0)
                 .lec.DeleteStream(paste0("mibrrStream", c(0 : nTargets)))
                 
                 ## Re-set the user's stream as the active RNG:
                 check <- length(userRng) > 0 & userRng != ""
                 if(check) .lec.CurrentStream(userRng)
             },
             
################################# ACCESSORS ####################################
             
             dataNames    = function() { colnames(data)                       },
             targets      = function() { targetVars                           },
             countMissing = function() { missCounts                           },
             
###--------------------------------------------------------------------------###
             
             getControl = function () {
                 "Return the 'control' list"
                 list(convThresh        = convThresh,
                      usePcStarts       = usePcStarts,
                                        #scale             = scale,
                      minPredCor        = minPredCor,
                      miceIters         = miceIters,
                      miceRidge         = miceRidge,
                      miceMethod        = miceMethod,
                      preserveStructure = preserveStructure,
                      checkConv         = checkConv,
                      optTraceLevel     = optTraceLevel,
                      optCheckKkt       = optCheckKkt,
                      optMethod         = optMethod,
                      optBoundLambda    = optBoundLambda)
             },
             
###--------------------------------------------------------------------------###
             
             getImpDataset = function(reset = FALSE) {
                 "Fill missing values to produce a single imputed dataset"
                 if(reset)
                     impRowsPool <<- 1 : sampleSizes[[length(sampleSizes)]][2]
                 
                 tmp <- data # Make a local copy of 'data'
                 
                 ## Randomly choose a posterior draw to use as imputations:
                 impRow      <-  sample(impRowsPool, 1)
                 impRowsPool <<- setdiff(impRowsPool, impRow)

                 for(j in targetVars) {
                     impSam                <- gibbsOut[[j]]$imps[impRow, ]
                     tmp[missList[[j]], j] <- impSam
                 }
                 
                 ## Restructure imputed data to match the raw data layout:
                 if(preserveStructure & length(ignoreVars) > 0)
                     data.frame(tmp, ignoredColumns)[ , rawNames]
                 else if(preserveStructure)
                     tmp[ , rawNames]
                 else
                     tmp
             },
          
######################### COMPLEX METHODS/SUBROUTINES ##########################
             
             resetMissing = function() {
                 "Replace imputed values with the original missing data code"
                 for(v in colnames(data)) {
                     ## Rebase indices as per R's preference :
                     missList[[v]]          <<- missList[[v]] + 1
                     data[missList[[v]], v] <<- missCode
                 }
             },
             
###--------------------------------------------------------------------------###

             processInputs = function() {
                 "Process and check the user inputs and isolate a set of target variables"
                 
                 ## Process the Raw Inputs ##
                 
                 ## Save the original variable names:
                 rawNames <<- colnames(data)
               
                 ## Set aside the 'ignored' columns:
                 if(length(ignoreVars) > 0) {
                     ignoredColumns <<- as.data.frame(data[ , ignoreVars])
                     if(length(ignoreVars) == 1)
                         colnames(ignoredColumns) <<- ignoreVars
                     data <<- data[ , setdiff(colnames(data), ignoreVars)]
                 }
                                 
                 ## Store some useful metadata:
                 nVar   <<- as.integer(ncol(data))
                 nPreds <<- as.integer(nVar - 1)
                 nObs   <<- as.integer(nrow(data))

                 if(doMcem) {
                     smoothingWindow <<- as.integer(
                         min(10, ceiling(iterations[1] / 10))
                     )
                     
                     totalIters <<- sum(iterations)
                 }
                 
                 rHats <<- list()
                 for(j in targetVars) rHats[[j]] <<- list()
                                                         
                 ## Replace missCode entries in data with NAs
                 if(!is.na(missCode)) {
                     userMissCode           <<- TRUE
                     data[data == missCode] <<- NA
                 }
                 else 
                     userMissCode <<- FALSE
                 
                 ## Generate a pool of potential imputation indices:
                 if(doImp)
                     impRowsPool <<- 1 : sampleSizes[[length(sampleSizes)]][2]

                 ## Check Inputs' Admissability ##
                 
                 ## Check for target variables. When no targets are given, all
                 ## incomplete variables not listed in 'ignoreVars' are imputed.
                 if(length(targetVars) == 0) {
                     if(doImp) {
                         targetCandidates <-
                             colnames(data)[!colnames(data) %in% ignoreVars]
                         warning("You did not specify any target variables, so I will impute the missing data on\nevery variable in 'data' that is not listed in 'ignoreVars'.\n",
                                 call.      = FALSE,
                                 immediate. = TRUE)        
                     } else {
                         stop("Please specify a DV.")
                     }
                 } else {
                     targetCandidates <- targetVars
                 }
                 
                 ## Make sure 'data' contains missing data that we can find:
                 rMat <- is.na(data)
                 if(userMissCode)                      
                     if(!any(rMat, na.rm = TRUE))
                         stop(paste0("The value you provided for 'missCode' (i.e., ",
                                     missCode,
                                     ") does not appear anywhere in 'data'.\nAre you sure that ",
                                     missCode,
                                     " encodes your missing data?\n")
                              )
                 
                 if(length(targetCandidates) > 1) 
                     completeTargets <- colMeans(rMat[ , targetCandidates]) == 0
                 else 
                     completeTargets <- mean(rMat[ , targetCandidates]) == 0
                 
                 if(doImp & all(completeTargets)) 
                     stop("Your target variables appear to be fully observed. Did you forget to provide a\nvalue for 'missCode'?\n")
                 
                 ## Select the final set of target variables:
                 if(doImp) {
                     targetVars <<- targetCandidates[!completeTargets]
                     if(any(completeTargets))
                         warning(paste0("The potential target variables {",
                                        paste(targetCandidates[completeTargets],
                                              collapse = ", "),
                                        "} are fully observed.\nThese items will not be imputed.\n"),
                                 call.      = FALSE,
                                 immediate. = TRUE)
                 }

                 ## Do some housekeeping that needs veridical 'targetVars' ##
                 
                 ## Re-order non-ignored data columns and store as data:
                 data <<- data.frame(
                     data[ , targetVars],
                     data[ , setdiff(colnames(data), targetVars)]
                 )
                 
                 ## Hack to deal with 1D matrix conversion to vector:
                 if(length(targetVars) == 1) colnames(data)[1] <<- targetVars
                                
                 ## Store nonresponse counts:
                 missCounts <<- sapply(colSums(is.na(data)), as.integer)

                 ## Create a list of missing elements in each target variable
                 ## NOTE: Subtract 1 from each index vector to base indices
                 ##       at 0 for C++
                 missList <<- lapply(data, function(x) which(is.na(x)) - 1)
                 
                 ## How many targets?
                 nTargets <<- as.integer(length(targetVars))

                 ## Make some starting value containers:
                 betaStarts <<- matrix(NA, nPreds + 1, nTargets)
                 tauStarts  <<- matrix(NA,  nPreds, nTargets)
                 
                 ## Initialize penalty parameter-related stuff:
                 lambda1Starts <<- rep(0.5, nTargets)
                 lambda2Starts <<- rep(nPreds / 10, nTargets)
                                 
                 if(doMcem) {
                     lambdaHistory <<-
                         lapply(targetVars,
                                function(x) {
                                    tmp <-
                                        matrix(NA, totalIters - 1, 2)
                                    colnames(tmp) <- c("lambda1", "lambda2")
                                    tmp
                                })
                     names(lambdaHistory) <<- targetVars
                 }
             },
             
###--------------------------------------------------------------------------###
             
             nameOutput = function() {
                 "Give the Gibb's sampling output pretty names"
                 if(checkConv) names(rHats) <<- targetVars
                 
                 for(v in targetVars) {
                     tmp <- setdiff(colnames(data), v)
                     
                     colnames(gibbsOut[[v]]$beta) <<- c("intercept", tmp)
                     colnames(gibbsOut[[v]]$tau)  <<- tmp
                 }
             },

###--------------------------------------------------------------------------###

             calcRHat = function(sims) {
                 "Compute a single split-chain Potential Scale Reduction Factor"
                 subChainLen <- floor(length(sims) / 2)
                 nSubChains  <- nChains * 2
                 
                 if(length(sims) %% nSubChains == 0) {
                     simsMat <- matrix(sims, ncol = nSubChains)
                 } else {
                     simsMat <- matrix(
                         sims[1 : (length(sims) - (nSubChains - 1))],
                         ncol = nSubChains
                     )
                 }
                 
                 wMean     <- colMeans(simsMat)
                 grandMean <- mean(simsMat)

                 wVar <- mean(apply(simsMat, 2, var))
                 bVar <- (subChainLen / (nSubChains - 1)) *
                     sum((wMean - grandMean)^2)
                 tVar <- ((subChainLen - 1) / subChainLen) * wVar +
                     (1 / subChainLen) * bVar
                 
                 sqrt(tVar / wVar)
             },

###--------------------------------------------------------------------------###
             
             computeRHats = function() {
                 "Compute the potential scale reduction factors"
                 for(j in targetVars) {
                     rHats[[j]]$beta  <<- apply(gibbsOut[[j]]$beta, 2, calcRHat)
                     rHats[[j]]$sigma <<- calcRHat(gibbsOut[[j]]$sigma)
                     
                     if(penalty != 0) {# Using a shrinkage prior?
                         rHats[[j]]$tau <<-
                             apply(gibbsOut[[j]]$tau, 2, calcRHat)
                         
                         if(!doMcem) {# Fully Bayesian estimation of lambdas?
                             if(penalty == 2)
                                 rHats[[j]]$lambda <<-
                                     apply(gibbsOut[[j]]$lambda, 2, calcRHat)
                             else
                                 rHats[[j]]$lambda <<-
                                     calcRHat(gibbsOut[[j]]$lambda[ , 1])
                         }# CLOSE if(!doMcem)
                     }# CLOSE if(penalty != 0)
                 }# END for(j in targetVars)
             },
             
###--------------------------------------------------------------------------###
             
             checkGibbsConv = function() {
                 "Check that the Gibb's sampler has converged everywhere"
                 for(j in targetVars) {
                     ## Find nonconvergent Gibbs samples:

#############################################################################################3
                     print(c("vanilla", "mibl", "mibrr")[penalty + 1])
                     print(" ")
                     print(j)
                     print(" ")
                     print("R-hats")
                     print(" ")
                     print(rHats[[j]]$beta)
                     print(" ")
                     print("Beta")
                     print(" ")
                     print(gibbsOut[[j]]$beta)
                     print(" ")
##############################################################################################
                     
                     badBetaCount <- sum(rHats[[j]]$beta > convThresh)
                     maxBetaRHat  <- max(rHats[[j]]$beta)
                     badSigmaFlag <- rHats[[j]]$sigma > convThresh
                     
                     if(penalty != 0) {
                         badTauCount  <- sum(rHats[[j]]$tau > convThresh)
                         maxTauRHat   <- max(rHats[[j]]$tau)

                         if(!doMcem) {
                             badLamCount <- sum(rHats[[j]]$lambda > convThresh)
                             maxLamRHat  <- max(rHats[[j]]$lambda)
                         }
                         else
                             badLamCount <- 0
                     }
                     else {
                         badTauCount <- 0
                         badLamCount <- 0
                     }
                     
                     ## Return warnings about possible failures of convergence:
                     if(badBetaCount > 0) {
                         warning(paste0("While imputing ",
                                        j,
                                        ", Beta's final Gibbs sample may not have converged.\n",
                                        badBetaCount,
                                        " R-Hats > ",
                                        convThresh,
                                        " with maximum R-Hat = ",
                                        round(maxBetaRHat, 4),
                                        ".\nConsider increasing the size of the (retained) Gibbs samples."),
                                 call.      = FALSE,
                                 immediate. = TRUE)
                     }
                     if(badTauCount > 0) {
                         warning(paste0("While imputing ",
                                        j,
                                        ", Tau's final Gibbs sample may not have converged.\n",
                                        badTauCount,
                                        " R-Hats > ",
                                        convThresh,
                                        " with maximum R-Hat = ",
                                        round(maxTauRHat, 4),
                                        ".\nConsider increasing the size of the (retained) Gibbs samples."),
                                 call.      = FALSE,
                                 immediate. = TRUE)
                     }
                     if(badSigmaFlag) {
                         warning(paste0("While imputing ",
                                        j,
                                        ", Sigma's final Gibbs sample ",
                                        "may not have converged.\nR-Hat = ",
                                        round(rHats[[j]]$sigma, 4),
                                        ".\nConsider increasing the size of the (retained) Gibbs samples."),
                                 call.      = FALSE,
                                 immediate. = TRUE)
                     }
                     if(badLamCount > 0) {
                         warning(paste0("While imputing ",
                                        j,
                                        ", Lambda's final Gibbs sample may not have converged.\n",
                                        badLamCount,
                                        " R-Hats > ",
                                        convThresh,
                                        " with maximum R-Hat = ",
                                        round(maxLamRHat, 4),
                                        ".\nConsider increasing the size of the (retained) Gibbs samples."),
                                 call.      = FALSE,
                                 immediate. = TRUE)
                     }
                 }# CLOSE for(j in targetVars)
             },
             
###--------------------------------------------------------------------------###

             startParams = function(restart = FALSE) {    
                 "Provide starting values for all parameters"

                 if(restart) {
                     for(j in 1 : nTargets) {
                         sigmaStarts[j]   <<- numMode(gibbsOut[[j]]$sigma)
                         tauStarts[ , j]  <<-
                             apply(gibbsOut[[j]]$tau, 2, numMode)
                         betaStarts[ , j] <<-
                             apply(gibbsOut[[j]]$beta, 2, numMode)
                     }
                     return()
                 }
                 
                 ## NOTE: We don't need real starting values for the intercepts.
                 ##       Their initial values will be sampled in the first
                 ##       iteration of the Gibbs sampler.
                 
                 ## Populate the starting values for Lambda:
                 if(penalty == 2) {
                     lambdaMat <<- cbind(
                         matrix(lambda1Starts, nTargets, 1),
                         matrix(lambda2Starts, nTargets, 1)
                     )
                 }
                 else if(penalty == 1) {
                     if(usePcStarts) getLambdaStarts()
                     lambdaMat <<-
                         cbind(matrix(lambda1Starts, nTargets, 1), 0)
                 }
                 else
                     lambdaMat <<- matrix(NA, nTargets, 2)
                 
                 rownames(lambdaMat) <<- targetVars
                 colnames(lambdaMat) <<- c("lambda1", "lambda2")
                 
                 ## Populate starting values for betas, taus, and sigma:
                 sigmaStarts <<- apply(data, 2, sd)[targetVars]
                 
                 for(j in 1 : nTargets) {
                     if(penalty == 2) {     # We're doing BEN
                         lam1 <- lambdaMat[j, 1]
                         lam2 <- lambdaMat[j, 2]
                         
                         tauPriorScale <- (8 * lam2 * sigmaStarts[j]) / lam1^2
                      
                         for(k in 1 : nPreds) {
                             tauDraw <- 0.0
                             while(tauDraw < 1.0)
                                 tauDraw <- rgamma(n     = 1,
                                                   shape = 0.5,
                                                   scale = tauPriorScale)
                             tauStarts[k, j] <<- tauDraw
                         }
                         
                         betaPriorCov <- diag(
                             1 / ((lam2 / sigmaStarts[j]) *
                                  (tauStarts[ , j] / (tauStarts[ , j] - 1.0))
                             )
                         )
                     }
                     else if(penalty == 1) {# We're doing BL
                         lam <- lambdaMat[j, 1]
                         
                         tauStarts[ , j] <<- rexp(nPreds, rate = (0.5 * lam^2))
                         
                         betaPriorCov <- sigmaStarts[j] * diag(tauStarts[ , j])
                     }
                     else                   # We're doing basic ridge 
                         betaPriorCov <- diag(rep(sigmaStarts[j], nPreds)) 
                     
                     betaStarts[ , j] <<-
                         c(0, rmvnorm(1, rep(0, nPreds), betaPriorCov))
                 }# CLOSE for(j in 1 : nTargets)
             },
             
###--------------------------------------------------------------------------###
             
             simpleImpute = function(intern = TRUE) { 
                 "Initially fill the missing values via single imputation"
                 cn     <- dataNames()
                 rFlags <- (missCounts > 0)[cn]
                 
                 ## Don't try to impute fully observed targets:
                 if(!any(rFlags)) return()
                 
                 ## Construct a predictor matrix for mice() to use:
                 predMat <-
                     quickpred(data, mincor = minPredCor, include = targetVars)
                 
                 ## Construct a vector of elementary imputation methods:
                 methVec         <- rep("", ncol(data))
                 methVec[rFlags] <- miceMethod
          
                 ## Singly impute the missing values:
                 miceOut <- mice(data            = data,
                                 m               = 1,
                                 maxit           = miceIters,
                                 method          = methVec,
                                 predictorMatrix = predMat,
                                 printFlag       = FALSE,
                                 ridge           = miceRidge)
                 
                 ## Replace missing values with their imputations:
                 if(intern)
                     setData(dat1 = mice::complete(miceOut, 1),
                             vars = cn[rFlags])
                 else
                     mice::complete(miceOut, 1) 
             },
             
###--------------------------------------------------------------------------###
             
             getLambdaStarts = function(nSamples = 25) {
                 "Use a variant of the method recommended by Park and Casella (2008) to get starting values for the MIBL penalty parameters"
                 
                 ## Fill any missing data with rough guesses:
                 impData <- simpleImpute(intern = FALSE)
                 
                 for(i in 1 : nTargets) {
                     check <- 0.90 * nObs > nPreds 
                     if(check) {# P << N
                         lmOut            <-  lm(impData[ , i] ~ impData[ , -i])
                         lambda1Starts[i] <<- nPreds * sd(resid(lmOut)) /
                             sum(abs(coef(lmOut)[-1]))
                     } else {
                         ## If P ~ N or  P > N, subsample data's columns and
                         ## repeatedly apply the Park & Casella (2008) method.
                         tmpLambda <- rep(NA, nSamples)
                         predCount <- round(0.90 * nObs)
                         for(j in 1 : nSamples) {
                             predSelector <-
                                 sample(c(1 : nVar)[-i], size = predCount)
                             tmpData      <- impData[ , predSelector]
                             lmOut        <- lm(impData[ , i] ~ tmpData)
                             tmpLambda[j] <- predCount * sd(resid(lmOut)) / 
                                 sum(abs(coef(lmOut)[-1]))
                         }# END for(j in 1 : nSamples)
                         lambda1Starts[i] <<- mean(tmpLambda)
                     }    
                 }
             }             
             
         )# END MibrrFit$methods()
