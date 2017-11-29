### Title:    MibrrFit Reference Class Definition
### Author:   Kyle M. Lang
### Created:  2017-NOV-28
### Modified: 2017-NOV-29
### Note:     MibrrFit is the metadata class for the MIBRR package

##--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------##
##  Copyright (C) 2017 Kyle M. Lang <k.m.lang@uvt.nl>                          ##  
##                                                                             ##
##  This file is part of MIBRR.                                                ##
##                                                                             ##
##  This program is free software: you can redistribute it and/or modify it    ##
##  under the terms of the GNU General Public License as published by the      ##
##  Free Software Foundation, either version 3 of the License, or (at you      ##
##  option) any later version.                                                 ##
##                                                                             ##
##  This program is distributed in the hope that it will be useful, but        ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of                 ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General   ##
##  Public License for more details.                                           ##
##                                                                             ##
##  You should have received a copy of the GNU General Public License along    ##
##  with this program. If not, see <http://www.gnu.org/licenses/>.             ##
##-----------------------------------------------------------------------------##

MibrrFit <- setRefClass("MibrrFit",
                        fields = list(
                            data              = "data.frame",
                            targetVars        = "character",
                            ignoreVars        = "character",
                            nImps             = "integer",
                            iterations        = "integer",
                            sampleSizes       = "list",
                            missCode          = "integer",
                            seed              = "integer",
                            doImp             = "logical",
                            doMcem            = "logical",
                            doBl              = "logical",
                            returnConvInfo    = "logical",
                            returnParams      = "logical",
                            verbose           = "logical",
                            convThresh        = "numeric",
                            lambda1Starts     = "numeric",
                            lambda2Starts     = "numeric",
                            usePcStarts       = "logical",
                            smoothingWindow   = "integer",
                            center            = "logical",
                            scale             = "logical",
                            adaptScales       = "logical",
                            simpleIntercept   = "logical",
                            minPredCor        = "numeric",
                            miceIters         = "integer",
                            miceRidge         = "numeric",
                            miceMethod        = "character",
                            fimlStarts        = "logical",
                            preserveStructure = "logical",
                            optTraceLevel     = "integer",
                            optCheckKkt       = "logical",
                            optMethod         = "character",
                            optBoundLambda    = "logical",
                            dataMeans         = "numeric",
                            dataScales        = "numeric",
                            gibbsOut          = "list",
                            ignoredColumns    = "data.frame",
                            rawNames          = "character",
                            impRowsPool       = "integer",
                            missList          = "list",
                            nChains           = "integer",
                            rHats             = "list",
                            lambdaMat         = "matrix",
                            betaStarts        = "matrix",
                            tauStarts         = "matrix",
                            sigmaStarts       = "numeric",
                            userMissCode      = "logical",
                            missCounts        = "integer",
                            nTargets          = "integer",
                            nVars             = "integer",
                            nPreds            = "integer",
                            nObs              = "integer"
                        )
                        )


MibrrFit$methods(

################################ CONSTRUCTOR ####################################
             
             initialize =
                 function(data              = data.frame(NULL),
                          targetVars        = "",
                          ignoreVars        = "",
                          nImps             = as.integer(NA),
                          iterations        = as.integer(NA),
                          sampleSizes       = list(),
                          missCode          = as.integer(NA),
                          seed              = as.integer(NA),
                          doImp             = as.logical(NA),
                          doMcem            = as.logical(NA),
                          doBl              = as.logical(NA),
                          returnConvInfo    = as.logical(NA),
                          returnParams      = as.logical(NA),
                          verbose           = as.logical(NA),
                          convThresh        = 1.1,
                          lambda1Starts     = as.numeric(NA),
                          lambda2Starts     = as.numeric(NA),
                          usePcStarts       = FALSE,
                          smoothingWindow   = as.numeric(NA),
                          center            = TRUE,
                          scale             = TRUE,
                          adaptScales       = TRUE,
                          simpleIntercept   = FALSE,
                          minPredCor        = 0.3,
                          miceIters         = 10L,
                          miceRidge         = 1e-4,
                          miceMethod        = "pmm",
                          fimlStarts        = FALSE,
                          preserveStructure = TRUE,
                          optTraceLevel     = 0L,
                          optCheckKkt       = TRUE,
                          optMethod         = "L-BFGS-B",
                          optBoundLambda    = TRUE,
                          dataMeans         = as.numeric(NA),
                          dataScales        = as.numeric(NA),
                          gibbsOut          = list(),
                          ignoredColumns    = data.frame(NULL),
                          rawNames          = "",
                          impRowsPool       = as.integer(NA),
                          missList          = list(),
                          nChains           = 1L,
                          rHats             = list(),
                          lambdaMat         = matrix(NA),
                          betaStarts        = matrix(NA),
                          tauStarts         = matrix(NA),
                          sigmaStarts       = as.numeric(NA),
                          userMissCode      = as.logical(NA),
                          missCounts        = as.integer(NA),
                          nTargets          = as.integer(NA),
                          nVars             = as.integer(NA),
                          nPreds            = as.integer(NA),
                          nObs              = as.integer(NA)
                          )
                 {
                     "Initialize an object of class MibrrFit"
                     data              <<- data
                     targetVars        <<- targetVars
                     ignoreVars        <<- ignoreVars
                     nImps             <<- nImps
                     iterations        <<- iterations
                     sampleSizes       <<- sampleSizes
                     missCode          <<- missCode
                     seed              <<- seed
                     doImp             <<- doImp
                     doMcem            <<- doMcem
                     doBl              <<- doBl
                     returnConvInfo    <<- returnConvInfo
                     returnParams      <<- returnParams
                     verbose           <<- verbose
                     convThresh        <<- convThresh
                     usePcStarts       <<- usePcStarts
                     center            <<- center
                     scale             <<- scale
                     adaptScales       <<- adaptScales
                     simpleIntercept   <<- simpleIntercept 
                     minPredCor        <<- minPredCor
                     miceIters         <<- miceIters
                     miceRidge         <<- miceRidge
                     miceMethod        <<- miceMethod
                     fimlStarts        <<- fimlStarts
                     preserveStructure <<- preserveStructure
                     optTraceLevel     <<- optTraceLevel
                     optCheckKkt       <<- optCheckKkt
                     optMethod         <<- optMethod
                     optBoundLambda    <<- optBoundLambda
                     dataMeans         <<- dataMeans
                     dataScales        <<- dataScales
                     gibbsOut          <<- gibbsOut
                     ignoredColumns    <<- ignoredColumns
                     rawNames          <<- rawNames
                     impRowsPool       <<- impRowsPool
                     missList          <<- missList
                     nChains           <<- nChains
                     lambdaMat         <<- lambdaMat
                     betaStarts        <<- betaStarts
                     tauStarts         <<- tauStarts
                     sigmaStarts       <<- sigmaStarts
                     
                     ## Save the original variable names:
                     rawNames <<- colnames(data)
                     
                     ## Set aside the 'ignored' columns:
                     ignoredColumns <<- as.data.frame(data[ , ignoreVars])
                     if(length(ignoreVars) == 1)
                         colnames(ignoredColumns) <<- ignoreVars
                   
                     ## Re-order non-ignored data columns and store as data
                     data <<- data.frame(
                         data[ , targetVars],
                         data[ , setdiff(colnames(data),
                                         c(targetVars, ignoreVars)
                                         )
                              ]
                     )
                     
                     ## Hack to deal with 1D matrix conversion to vector:
                     if(length(targetVars) == 1)
                         colnames(data)[1] <<- targetVars
                     
                     ## Store some useful metadata:
                     nTargets <<- nt <- as.integer(length(targetVars))
                     nVars    <<- nv <- as.integer(
                                      ncol(data) - length(ignoreVars)
                                  )
                     nPreds   <<- np <- as.integer(nv - 1)
                     nObs     <<- as.integer(nrow(data))
                    
                     tauStarts <<- betaStarts <<- matrix(NA, np, nt)
                                         
                     lambda1Starts <<- rep(0.5, nt)
                     lambda2Starts <<- rep(np / 10, nt)
                                         
                     smoothingWindow <<- as.integer(
                         min(10, ceiling(iterations[1] / 10))
                     )
                     
                     rHats <<- list()
                     for(j in targetVars)
                         rHats[[j]] <<- list(beta = NA, tau = NA, sigma = NA)
                   
                     ## Replace missCode entries in data with NAs
                     if(!is.na(missCode)) {
                         userMissCode           <<- TRUE
                         data[data == missCode] <<- NA
                     }
                     else {
                         userMissCode <<- FALSE
                     }

                     missCounts <<- sapply(colSums(is.na(data)), as.integer)
                     
                     ## Create a list of missing elements in each target variable
                     ## NOTE: Subtract 1 from each index vector to base indices
                     ##       at 0 for C++
                     missList <<- lapply(data, function(x) which(is.na(x)) - 1)
                 },
             
################################### MUTATORS ####################################

             setDataMeans  = function(dataMeans)  { dataMeans  <<- dataMeans   },
             setDataScales = function(dataScales) { dataScales <<- dataScales  },

###---------------------------------------------------------------------------###
             
             setData = function(data) { data[ , colnames(data)] <<- data       },
             
###---------------------------------------------------------------------------###
             
             setControl = function(x) {
                 "Assign the control parameters"
                 ints <- c("smoothingWindow", "miceIters", "optTraceLevel")
                 
                 for(n in names(x)) {
                     if(n %in% ints) field(n, as.integer(x[[n]]))
                     else            field(n, x[[n]])
                 }
             },

################################# ACCESSORS #####################################

             dataNames    = function() { colnames(data)                        },
             targets      = function() { targetVars                            },
             countMissing = function() { missCounts                            },
             
###---------------------------------------------------------------------------###
             
             getControl = function () {
                 "Return the 'control' list"
                 list(convThresh        = convThresh,
                      usePcStarts       = usePcStarts,
                      center            = center,
                      scale             = scale,
                      adaptScales       = adaptScales,
                      simpleIntercept   = simpleIntercept, 
                      minPredCor        = minPredCor,
                      miceIters         = miceIters,
                      miceRidge         = miceRidge,
                      miceMethod        = miceMethod,
                      fimlStarts        = fimlStarts,
                      preserveStructure = preserveStructure,
                      optTraceLevel     = optTraceLevel,
                      optCheckKkt       = optCheckKkt,
                      optMethod         = optMethod,
                      optBoundLambda    = optBoundLambda)
             },
             
###---------------------------------------------------------------------------###
             
             getImpDataset = function() {
                 "Fill missing values to produce a single imputed dataset"
                 tmp <- data # Make a local copy of 'data'
                 
                 ## Randomly choose a posterior draw to use as imputations:
                 impRow      <-  sample(impRowsPool, 1)
                 impRowsPool <<- setdiff(impRowsPool, impRow)
                 
                 for(j in targetVars) {
                     impSam <- gibbsOut[[j]]$imps[impRow, ]
                     tmp[missList[[j]], j] <- impSam + dataMeans[j]
                 }
                 
                 ## Restructure imputed data to match the raw data layout:
                 if(preserveStructure)
                     data.frame(tmp, ignoredColumns)[ , rawNames]
                 else
                     tmp
             },
          
########################## COMPLEX METHODS/SUBROUTINES ##########################
             
             applyMissCode = function() {
                 "Construct an integer-valued missing data code"
                 if(is.na(missCode)) {
                     if(max(abs(data), na.rm = TRUE) < 1.0) {
                         missCode <<- -9
                     }
                     else {
                         codeMag   <-
                             floor(log10(max(abs(data), na.rm = TRUE))) + 2
                         missCode <<- -(10^codeMag - 1)
                     }
                 }
                 data[is.na(data)] <<- missCode
             },

###---------------------------------------------------------------------------###
             
             checkInputs = function() {
                 "Check the user inputs and isolate a set of target variables"
                 
                 ## Check for target variables. When no targets are given, all
                 ## incomplete variables not listed in 'ignoreVars' are imputed.
                 if(is.null(targetVars)) {
                     if(doImp) {
                         targetCandidates <-
                             colnames(data)[!colnames(data) %in% ignoreVars]
                         warning("You did not specify any target variables, so I will impute the missing data on\nevery variable in 'data' that is not listed in 'ignoreVars'.\n")        
                     } else {
                         stop("Please specify a DV.")
                     }
                 } else {
                     targetCandidates <- targetVars
                 }
                 
                 ## Make sure 'data' contains missing data that we can find:
                 if(is.na(missCode)) {
                     rMat <- is.na(data)
                 } else {
                     rMat <- data == missCode
                     
                     if(!any(rMat, na.rm = TRUE))
                         stop(paste0("The value you provided for 'missCode' (i.e., ",
                                     missCode,
                                     ") does not appear anywhere in 'data'.\nAre you sure that ",
                                     missCode,
                                     " encodes your missing data?\n")
                              )
                 }
                 
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
                         warning(
                             paste0("The potential target variables {",
                                    paste(targetCandidates[completeTargets],
                                          collapse = ", "),
                                    "} are fully observed.\nThese items will not be imputed.\n")
                         )
                 }
             },

###---------------------------------------------------------------------------###
             
             scaleData = function(revert = FALSE, compStats = TRUE) {
                 "Standardize the columns of data"
                 if(!revert) {# Doing initial scaling
                     if(compStats) {# Compute summary stats
                         if(scale)
                             dataScales <<- unlist(lapply(data, sd))
                         else
                             dataScales <<- rep(1, nVar)
                         
                         ## Mean center data:
                         if(center) {
                             dataMeans <<- colMeans(data)
                             data      <<- as.data.frame(
                                 scale(data, center = TRUE, scale = FALSE)
                             )
                         } else {
                             dataMeans <<- rep(0, nVar)
                         }
                         
                         names(dataMeans) <<-
                             names(dataScales) <<- colnames(data)
                     }
                     else {# Don't re-compute summary stats
                         data <<- data -
                             data.frame(
                                 matrix(dataMeans, nObs, nVar, byrow = TRUE)
                             )
                     }
                 }
                 else {# Reverting the data to its original scaling
                     data <<-
                         data + data.frame(
                                    matrix(dataMeans, nObs, nVar, byrow = TRUE)
                                )
                 }
             },
             
###---------------------------------------------------------------------------###

             nameOutput = function() {
                 "Give the Gibb's sampling output pretty names"
                 if(returnConvInfo) names(rHatList) <- targetVars
                 
                 if(returnParams)
                     for(v in targetVars) {
                         tmp <- setdiff(colnames(data), v)
                         
                         colnames(gibbsOut[[v]]$beta) <<- c("intercept", tmp)
                         colnames(gibbsOut[[v]]$tau)  <<- tmp
                     }
             },


###---------------------------------------------------------------------------###

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

###---------------------------------------------------------------------------###
             
             computeRHats = function() {
                 "Compute the potential scale reduction factors"
                 for(j in targetVars) {
                     rHats[[j]]$beta  <<- apply(gibbsOut[[j]]$beta, 2, calcRHat)
                     rHats[[j]]$tau   <<- apply(gibbsOut[[j]]$tau,  2, calcRHat)
                     rHats[[j]]$sigma <<- calcRHat(gibbsOut[[j]]$sigma)
                 }
             },

###---------------------------------------------------------------------------###
             
             checkGibbsConv = function() {
                 "Check that the Gibb's sampler has converged everywhere"
                 for(j in targetVars) {
                     ## Find nonconvergent Gibbs samples:
                     badBetaCount <- sum(rHats[[j]]$beta > convThresh)
                     maxBetaRHat  <- max(rHats[[j]]$beta)
                     badTauCount  <- sum(rHats[[j]]$tau > convThresh)
                     maxTauRHat   <- max(rHats[[j]]$tau)
                     badSigmaFlag <- rHats[[j]]$sigma > convThresh
                     
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
                                        ".\nConsider increasing the size of the (retained) Gibbs samples.")
                                 )
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
                                        ".\nConsider increasing the size of the (retained) Gibbs samples.")
                                 )
                     }
                     if(badSigmaFlag) {
                         warning(paste0("While imputing ",
                                        j,
                                        ", Sigma's final Gibbs sample ",
                                        "may not have converged.\nR-Hat = ",
                                        round(sigmaRHat, 4),
                                        ".\nConsider increasing the size of the (retained) Gibbs samples."))
                     }
                 }# CLOSE for(j in targetVars)
             },
             
###---------------------------------------------------------------------------###

             startParams = function(restart = FALSE) {    
                 "Provide starting values for all parameters"
                 
                 if(restart) {
                     for(j in 1 : nTargets) {
                         sigmaStarts[j]   <<- mean(gibbsOut[[j]]$sigma)
                         tauStarts[ , j]  <<- colMeans(gibbsOut[[j]]$tau)
                         betaStarts[ , j] <<-
                             colMeans(gibbsOut[[j]]$beta[ , -1])
                     }
                     return()
                 }
                 
                 ## NOTE: We don't need to start the intercept. It's initial
                 ##       value will be sampled in the first iteration of the
                 ##       Gibbs sampler.
                 
                 ## Populate the starting values for Lambda:
                 if(!doBl)
                     lambdaMat <<- cbind(
                         matrix(lambda1Starts, nTargets, 1),
                         matrix(lambda2Starts, nTargets, 1)
                     )
                 else 
                     lambdaMat <<- cbind(matrix(lambda1Starts, nTargets, 1), 0)
                 
                 ## Populate starting values for betas, taus, and sigma:
                 sigmaStarts <<- dataScales[targetVars]
                                
                 for(j in 1 : nTargets) {
                     if(!doBl) {
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

                         betaStarts[ , j] <<-
                             rmvnorm(1, rep(0, nPreds), betaPriorCov)
                     }
                     else {# We're doing BL
                         lam <- lambdaMat[j, 1]
                         
                         tauStarts[ , j] <<- rexp(nPreds, rate = (0.5 * lam^2))
                         
                         betaPriorCov <- sigmaStarts[j] * diag(tauStarts[ , j])
                         betaStarts[ , j] <<-
                             rmvnorm(1, rep(0, nPreds), betaPriorCov)
                     }
                 }# CLOSE for(j in 1 : nTargets)
             },

###---------------------------------------------------------------------------###
             
             smoothLambda = function() {
                 i     <- iterations[1]
                 range <- (i - smoothingWindow + 1) : i

                 for(j in targetVars)
                     lambdaMat[j, ] <<- colMeans(lambdaHistory[[j]][range, ])
             }
             
         )# END MibrrFit$methods()
