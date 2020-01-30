### Title:    MibrrFit Reference Class Definition
### Author:   Kyle M. Lang
### Created:  2017-11-28
### Modified: 2020-01-30
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
                            data         = "data.frame",
                            targetVars   = "character",
                            ignoreVars   = "character",
                            iterations   = "integer",
                            sampleSizes  = "list",
                            missCode     = "integer",
                            seed         = "ANY",
                            doImp        = "logical",
                            doMcem       = "logical",
                            verbose      = "logical",
                                        #l1Pars       = "numeric",
                                        #l2Pars       = "numeric",
                            ignoredCols  = "data.frame",
                            rawNames     = "character",
                            impRowsPool  = "matrix",
                            missList     = "list",
                            nChains      = "integer",
                            rHats        = "list",
                            mct          = "list",
                            userMissCode = "logical",
                            missCounts   = "integer",
                            nTargets     = "integer",
                            nVar         = "integer",
                            nPreds       = "integer",
                            nObs         = "integer",
                                        #totalIters   = "integer",
                            rng0         = "character",
                            userRng      = "character",
                            streams      = "character",
                            ridge        = "numeric",
                            penalty      = "integer",
                            chains       = "list",
                            nCores       = "integer",
                            control      = "list"
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
                          penalty     = 2L,
                          nChains     = 1L,
                          nCores      = 1L)
                 {
                     "Initialize an object of class MibrrFit"
                     data        <<- data
                     targetVars  <<- targetVars
                     ignoreVars  <<- ignoreVars
                     iterations  <<- iterations
                     sampleSizes <<- sampleSizes
                     missCode    <<- missCode
                     doImp       <<- doImp
                     doMcem      <<- doMcem
                     verbose     <<- verbose
                     nChains     <<- nChains
                     seed        <<- seed
                     userRng     <<- userRng
                     ridge       <<- ridge
                     penalty     <<- penalty
                     nCores      <<- nCores
                 },
             
################################### MUTATORS ###################################
             
             setData = function(dat1, vars = NULL) {
                 if(is.null(vars))
                     data[ , colnames(dat1)] <<- dat1
                 else
                     data[ , vars] <<- dat1[ , vars]
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
                 streams <<- c("mibrrStream0",
                               paste0("c",
                                      rep(1 : nChains, each = nTargets),
                                      targetVars)
                               )
                 
                 .lec.CreateStream(streams)
                 
                 ## Set the 'master' RNG stream:
                 rng0 <<- .lec.CurrentStream("mibrrStream0")
             },

###--------------------------------------------------------------------------###

             cleanRng = function() {
                 "Clean up the RNG state"
                 .lec.CurrentStreamEnd(rng0)
                 .lec.DeleteStream(streams)
                 
                 ## Re-set the user's stream as the active RNG:
                 check <- length(userRng) > 0 & userRng != ""
                 if(check) .lec.CurrentStream(userRng)
             },

###--------------------------------------------------------------------------###
             
################################# ACCESSORS ####################################
             
             getImpDataset = function(reset = FALSE) {
                 "Fill missing values to produce a single imputed dataset"
                 if(reset)
                     impRowsPool <<- as.matrix(
                         expand.grid(chain = 1 : nChains,
                                     row   = 1 : sampleSizes[[length(sampleSizes)]][2]
                                     )
                     )
                 
                 tmp <- data # Make a local copy of 'data'
                 
                 ## Randomly choose a posterior draw to use as imputations:
                 index       <-  sample(1 : nrow(impRowsPool), 1)
                 impRow      <-  impRowsPool[index, ]
                 impRowsPool <<- impRowsPool[-index, ]
                 
                 for(j in targetVars) {
                     impSam <-
                         chains[[impRow["chain"]]]$parameters[[j]]$imps[impRow["row"], ]
                     tmp[missList[[j]], j] <- impSam
                 }
                 
                 ## Restructure imputed data to match the raw data layout:
                 if(length(ignoreVars) > 0)
                     data.frame(tmp, ignoredCols)[rawNames]
                 else
                     tmp[rawNames]
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
                     ignoredCols <<- data[ignoreVars]
                     data        <<- data[setdiff(colnames(data), ignoreVars)]
                 }
                 
                 ## Store some useful metadata:
                 nVar   <<- as.integer(ncol(data))
                 nPreds <<- as.integer(nVar - 1)
                 nObs   <<- as.integer(nrow(data))

                                        #if(doMcem) totalIters <<- sum(iterations)
                 
                                        #rHats <<- list() 
                                        #for(j in targetVars) rHats[[j]] <<- list()
                 
                 ## Replace missCode entries in data with NAs
                 if(!is.na(missCode)) {
                     userMissCode           <<- TRUE
                     data[data == missCode] <<- NA
                 }
                 else 
                     userMissCode <<- FALSE
                 
                 ## Generate a pool of potential imputation indices:
                 if(doImp)
                     impRowsPool <<- as.matrix(
                         expand.grid(chain = 1 : nChains,
                                     row   = 1 : sampleSizes[[length(sampleSizes)]][2]
                                     )
                     )
                 
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
                     data[targetVars],
                     data[setdiff(colnames(data), targetVars)]
                 )
                 
                 ## Store nonresponse counts:
                 missCounts <<- sapply(colSums(is.na(data)), as.integer)
                 
                 ## Create a list of missing elements in each target variable
                 ## NOTE: Subtract 1 from each index vector to base indices
                 ##       at 0 for C++
                 missList <<- lapply(data, function(x) which(is.na(x)) - 1)
                 
                 ## How many targets?
                 nTargets <<- as.integer(length(targetVars))

                 ## Prep containers for R-Hats and MCTs:
                 rHats <<- mct <<- vector("list", nTargets)
                 names(rHats) <<- names(mct) <<- targetVars
             },
             
###--------------------------------------------------------------------------###
             
             initChains = function() {
                 "Initialize the 'MibrrChain' objects"
                 for(k in 1 : nChains) {
                     chains[[k]] <<- MibrrChain(chain       = k,
                                                data        = data,
                                                targetVars  = targetVars,
                                                iterations  = iterations,
                                                sampleSizes = sampleSizes,
                                                missList    = missList,
                                                doMcem      = doMcem,
                                                verbose     = verbose,
                                                ridge       = ridge,
                                                penalty     = penalty,
                                                control     = control)
                     
                     ## Set control parameters for the 'MibrrChain' objects:
                     setControl(x = control, where = chains[[k]])
                 }
             },
             
###--------------------------------------------------------------------------###
             
             nameOutput = function() {
                 "Give the Gibb's sampling output pretty names"
                 if(control$checkConv) names(rHats) <<- targetVars

                 for(k in 1 : nChains)
                     for(j in targetVars) {
                         tmp <- setdiff(colnames(data), j)
                         
                         colnames(chains[[k]]$parameters[[j]]$beta) <<- c("intercept", tmp)
                         colnames(chains[[k]]$parameters[[j]]$tau)  <<- tmp
                     }
             },

###--------------------------------------------------------------------------###
             
             getSamples = function(what, target) {
                 "Extract (unpooled) posterior samples"
                 out <- list()
                 for(k in 1 : nChains)
                     out[[k]] <- chains[[k]]$parameters[[target]][[what]]
                 
                 out
             },
             
###--------------------------------------------------------------------------###
             
                                        #splitSamples = function(what, target)
                                        #    unlist(
                                        #        lapply(getSamples(what, target), split),
                                        #        recursive = FALSE
                                        #    ),
             
###--------------------------------------------------------------------------###

             poolSamples = function(what, target) {
                 "Pool posterior samples across chains"
                 tmp <- list()
                 for(k in 1 : nChains)
                     tmp[[k]] <- chains[[k]]$parameters[[target]][[what]]

                 if(what %in% c("sigma", "lambda1", "lambda2"))
                     unlist(tmp, use.names = FALSE)
                 else
                     do.call(rbind, tmp)
             },
             
###--------------------------------------------------------------------------###
             
             calcRHats = function(sims, trans = FALSE, mv = TRUE) {
                 ## Prepare the samples for coda:
                 sims <- prepSams(sims)

                 ## Calculate the R-Hats:
                 gelman.diag(sims, transform = trans, multivariate = mv)
             },
             
###--------------------------------------------------------------------------###
             
             computeRHats = function() {
                 "Compute the potential scale reduction factors"
                 for(j in targetVars) {
                     tmp <- calcRHats(getSamples("sigma", j), trans = TRUE)
                     rHats[[j]]$sigma <<- tmp
                     
                     tmp <- calcRHats(getSamples("beta", j))
                     rHats[[j]]$beta <<- tmp
                     
                     if(penalty != 0) {# Using a shrinkage prior?
                         tmp <- calcRHats(getSamples("tau", j), trans = TRUE)
                         rHats[[j]]$tau <<- tmp
                         
                         if(!doMcem) {# Fully Bayesian estimation of lambdas?
                             tmp <- calcRHats(getSamples("lambda1", j),
                                              trans = TRUE)
                             rHats[[j]]$lambda1 <<- tmp
                             
                             if(penalty == 2) {
                                 tmp <- calcRHats(getSamples("lambda2", j),
                                                  trans = TRUE)
                                 rHats[[j]]$lambda2 <<- tmp
                             }
                         }# CLOSE if(!doMcem)
                     }# CLOSE if(penalty != 0)
                 }# END for(j in targetVars)
             },
             
###--------------------------------------------------------------------------###
             
             doMct = function(what, target, n) {
                 ## Prepare the samples:
                 sams <- getSamples(what, target)
                 sams <- lapply(sams, tail, n = n)
                 sams <- unlist(lapply(sams, split), recursive = FALSE)

                 ## Conduct the MCT:
                 mct(sams)
             },

###--------------------------------------------------------------------------###

             lambdaMct = function() {
                 for(j in targetVars) {
                     mct[[j]]$lambda1 <<- doMct("lambda1", j, iterations[3])
                     
                     if(penalty == 2)
                         mct[[j]]$lambda2 <<- doMct("lambda2", j, iterations[3])
                     
                     if(doMcem)
                         mct[[j]]$logLik  <<- doMct("logLik", j, iterations[3])
                 }
             },
             
###--------------------------------------------------------------------------###
             
             checkGibbsConv = function() {
                 "Check that the Gibb's sampler has converged everywhere"
                 for(j in targetVars) {
                     ## Find nonconvergent Gibbs samples:
                     badBetaCount <- sum(rHats[[j]]$beta$psrf[ , 2] > control$convThresh)
                     maxBetaRHat  <- max(rHats[[j]]$beta$psrf[ , 1])
                     badSigmaFlag <- rHats[[j]]$sigma$psrf[ , 2] > control$convThresh
                     
                     badTauCount <- badL1Count <- badL2Count <- 0
                     if(penalty != 0) {
                         badTauCount <- sum(rHats[[j]]$tau$psrf[ , 2] > control$convThresh)
                         maxTauRHat  <- max(rHats[[j]]$tau$psrf[ , 1])
                         
                         if(!doMcem) {
                             badL1Count <-
                                 sum(rHats[[j]]$lambda1$psrf[ , 2] > control$convThresh)
                             maxL1RHat  <- max(rHats[[j]]$lambda1$psrf[ , 1])
                             
                             if(penalty == 2) {
                                 badL2Count <-
                                     sum(rHats[[j]]$lambda2$psrf[ , 2] > control$convThresh)
                                 maxL2RHat  <- max(rHats[[j]]$lambda2$psrf[ , 1])
                             }
                         }                            
                     }
                     
                     ## Return warnings about possible failures of convergence:
                     if(badBetaCount > 0) {
                         warning(paste0("While imputing ",
                                        j,
                                        ", Beta's final Gibbs sample may not have converged.\n",
                                        badBetaCount,
                                        " R-Hats > ",
                                        control$convThresh,
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
                                        control$convThresh,
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
                                        round(rHats[[j]]$sigma$psrf[ , 1], 4),
                                        ".\nConsider increasing the size of the (retained) Gibbs samples."),
                                 call.      = FALSE,
                                 immediate. = TRUE)
                     }
                     if(badL1Count > 0) {
                         warning(paste0("While imputing ",
                                        j,
                                        ", Lambda1's final Gibbs sample may not have converged.\n",
                                        badL1Count,
                                        " R-Hats > ",
                                        control$convThresh,
                                        " with maximum R-Hat = ",
                                        round(maxL1RHat, 4),
                                        ".\nConsider increasing the size of the (retained) Gibbs samples."),
                                 call.      = FALSE,
                                 immediate. = TRUE)
                     }
                     if(badL2Count > 0) {
                         warning(paste0("While imputing ",
                                        j,
                                        ", Lambda2's final Gibbs sample may not have converged.\n",
                                        badL2Count,
                                        " R-Hats > ",
                                        control$convThresh,
                                        " with maximum R-Hat = ",
                                        round(maxL2RHat, 4),
                                        ".\nConsider increasing the size of the (retained) Gibbs samples."),
                                 call.      = FALSE,
                                 immediate. = TRUE)
                     }
                 }# CLOSE for(j in targetVars)
             },
             
###--------------------------------------------------------------------------###
             
             simpleImpute = function(intern = TRUE) { 
                 "Initially fill the missing values via single imputation"
                 cn     <- colnames(data)
                 rFlags <- (missCounts > 0)[cn]
                 
                 ## Don't try to impute fully observed targets:
                 if(!any(rFlags)) return()
                 
                 ## Construct a predictor matrix for mice() to use:
                 predMat <- mice::quickpred(data,
                                            mincor  = control$minPredCor,
                                            include = targetVars)
                 
                 ## Construct a vector of elementary imputation methods:
                 methVec         <- rep("", ncol(data))
                 methVec[rFlags] <- control$miceMethod
                 
                 ## Singly impute the missing values:
                 miceOut <- with(control,
                                 mice(data            = data,
                                      m               = 1,
                                      maxit           = miceIters,
                                      method          = methVec,
                                      predictorMatrix = predMat,
                                      printFlag       = FALSE,
                                      ridge           = miceRidge)
                                 )
                 
                 ## Replace missing values with their imputations:
                 if(intern)
                     setData(dat1 = mice::complete(miceOut, 1),
                             vars = cn[rFlags])
                 else
                     mice::complete(miceOut, 1) 
             }
             
###--------------------------------------------------------------------------###
             
                                        #getLambdaStarts = function(nSamples = 25) {
                                        #    "Use a variant of the method recommended by Park and Casella (2008) to get starting values for the MIBL penalty parameters"
                                        #    
                                        #    ## Fill any missing data with rough guesses:
                                        #    impData <- simpleImpute(intern = FALSE)
                                        #    
                                        #    for(i in 1 : nTargets) {
                                        #        check <- 0.90 * nObs > nPreds 
                                        #        if(check) {# P << N
                                        #            lmOut            <-  lm(impData[ , i] ~ impData[ , -i])
                                        #            lambda1Starts[i] <<- nPreds * sd(resid(lmOut)) /
                                        #                sum(abs(coef(lmOut)[-1]))
                                        #        } else {
                                        #            ## If P ~ N or  P > N, subsample data's columns and
                                        #            ## repeatedly apply the Park & Casella (2008) method.
                                        #            tmpLambda <- rep(NA, nSamples)
                                        #            predCount <- round(0.90 * nObs)
                                        #            for(j in 1 : nSamples) {
                                        #                predSelector <-
                                        #                    sample(c(1 : nVar)[-i], size = predCount)
                                        #                tmpData      <- impData[ , predSelector]
                                        #                lmOut        <- lm(impData[ , i] ~ tmpData)
                                        #                tmpLambda[j] <- predCount * sd(resid(lmOut)) / 
                                        #                    sum(abs(coef(lmOut)[-1]))
                                        #            }# END for(j in 1 : nSamples)
                                        #            lambda1Starts[i] <<- mean(tmpLambda)
                                        #        }    
                                        #    }
                                        #}             
             
         )# END MibrrFit$methods()
