### Title:    MibrrFit Reference Class Definition
### Author:   Kyle M. Lang
### Created:  2017-NOV-28
### Modified: 2017-NOV-28
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

MibrrFit <- setRefClass("MibrrFit")

MibrrFit$fields(
             _data              = "data.frame",
             _targetVars        = "character",
             _ignoreVars        = "character"
             _nImps             = "integer",
             _iterations        = "integer",
             _sampleSizes       = "list",
             _missCode          = "integer",
             _seed              = "integer",
             _doImp             = "logical",
             _doFullBayes       = "logical",
             _returnConvInfo    = "logical",
             _returnParams      = "logical",
             _verbose           = "logical",
             _convThresh        = "numeric",
             _lambda1Starts     = "numeric",
             _lambda2Starts     = "numeric",
             _usePcStarts       = "logical",
             _smoothingWindow   = "integer",
             _center            = "logical",
             _scale             = "logical",
             _adaptScales       = "logical",
             _simpleIntercept   = "logical",
             _minPredCor        = "numeric",
             _miceIters         = "integer",
             _miceRidge         = "numeric",
             _miceMethod        = "character",
             _fimlStarts        = "logical",
             _preserveStructure = "logical",
             _optTraceLevel     = "integer",
             _optCheckKkt       = "logical",
             _optMethod         = "character",
             _optBoundLambda    = "logical",
             _dataMeans         = "numeric",
             _dataScales        = "numeric",
             _gibbsOut          = "list",
             _ignoredColumns    = "data.frame",
             _rawNames          = "character",
             _impRowsPool       = "integer",
             _missList          = "list",
             _nChains           = "integer"
             _rHats             = "list",
             _lambdaMat         = "matrix",
             _betaStarts        = "matrix",
             _tauStarts         = "matrix",
             _sigmaStarts       = "numeric"
             _userMissCode      = "logical"
             _missCounts        = "integer"
         )# END MibrrFit$fields()


MibrrFit$methods(

################################ CONSTRUCTOR ####################################
             
             initialize =
                 function(data              = data.frame(NULL),
                          targetVars        = "",
                          ignoreVars        = ""
                          nImps             = as.integer(NA),
                          iterations        = as.integer(NA),
                          sampleSizes       = list(),
                          missCode          = as.integer(NA),
                          seed              = as.integer(NA),
                          doImp             = as.logical(NA),
                          doFullBayes       = as.logical(NA),
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
                          missCounts        = as.integer(NA)
                          )
                 {
                     "Initialize an object of class MibrrFit"
                     _data              <<- data
                     _targetVars        <<- targetVars
                     _ignoreVars        <<- ignoreVars
                     _nImps             <<- nImps
                     _iterations        <<- iterations
                     _sampleSizes       <<- sampleSizes
                     _missCode          <<- missCode
                     _seed              <<- seed
                     _doImp             <<- doImp
                     _doFullBayes       <<- doFullBayes
                     _returnConvInfo    <<- returnConvInfo
                     _returnParams      <<- returnParams
                     _verbose           <<- verbose
                     _convThresh        <<- convThresh
                     _usePcStarts       <<- usePcStarts
                     _center            <<- center
                     _scale             <<- scale
                     _adaptScales       <<- adaptScales
                     _simpleIntercept   <<- simpleIntercept 
                     _minPredCor        <<- minPredCor
                     _miceIters         <<- miceIters
                     _miceRidge         <<- miceRidge
                     _miceMethod        <<- miceMethod
                     _fimlStarts        <<- fimlStarts
                     _preserveStructure <<- preserveStructure
                     _optTraceLevel     <<- optTraceLevel
                     _optCheckKkt       <<- optCheckKkt
                     _optMethod         <<- optMethod
                     _optBoundLambda    <<- optBoundLambda
                     _dataMeans         <<- dataMeans
                     _dataScales        <<- dataScales
                     _gibbsOut          <<- gibbsOut
                     _ignoredColumns    <<- ignoredColumns
                     _rawNames          <<- rawNames
                     _impRowsPool       <<- impRowsPool
                     _missList          <<- missList
                     _nChains           <<- nChains
                     _lambdaMat         <<- lambdaMat
                     _betaStarts        <<- betaStarts
                     _tauStarts         <<- tauStarts
                     _sigmaStarts       <<- sigmaStarts
                     
                     ## Save the original variable names:
                     rawNames <<- colnames(data)
                     
                     ## Set aside the 'ignored' columns:
                     _ignoredColumns <<- as.data.frame(data[ , ignoreVars])
                     if(length(ignoreVars) == 1)
                         colnames(_ignoredColumns) <<- ignoreVars
                     
                     ## Re-order non-ignored data columns and store as _data
                     _data <<- data.frame(
                         data[ , targetVars],
                         data[ , setdiff(colnames(data),
                                         c(targetVars, ignoreVars)
                                         )
                              ]
                     )
                     
                     ## Hack to deal with 1D matrix conversion to vector:
                     if(length(targetVars) == 1)
                         colnames(_data)[1] <<- targetVars
                     
                     nTargets <- length(_targetVars)
                     nPreds   <- ncol(_data) - 1
                     
                     _lambda1Starts <<- rep(0.5, nTargets)
                     _lambda2Starts <<- rep(nPreds / 10, nTargets)

                     _smoothingWindow <<- min(10, ceiling(_iterations[1] / 10))

                     _rHats <<- list()
                     for(j in _targetVars)
                         _rHats[[j]] <- list(beta = NA, tau = NA, sigma = NA)

                     ## Replace missCode entries in _data with NAs
                     if(!is.na(missCode)) {
                         _userMissCode            <<- TRUE
                         _data[_data == missCode] <<- NA
                     } else {
                         _userMissCode <- FALSE
                     }
                     
                     ## Create a list of missing elements in each target variable
                     ## NOTE: Subtract 1 from each index vector to base indices
                     ##       at 0 for C++
                     _missList <<- lapply(_data, function(x) which(is.na(x)) - 1)

                     _missCounts <<- colSums(is.na(_data))
                 },

################################### MUTATORS ####################################

             setDataMeans  = function(dataMeans)  { _dataMeans  <<- dataMeans  },
             setDataScales = function(dataScales) { _dataScales <<- dataScales },

###---------------------------------------------------------------------------###
             
             setData = function(data) { _data[ , colnames(data)] <<- data      },
             
###---------------------------------------------------------------------------###
             
             setControl = function(x) {
                 "Assign the control parameters"
                 x    <- lapply(x, function(y) paste0("_", y))
                 ints <- c("_smoothingWindow", "_miceIters", "_optTraceLevel")
                 
                 for(n in names(x)) {
                     if(n %in% ints) field(n, as.integer(x[[n]]))
                     else            field(n, x[[n]])
                 }
             },

################################# ACCESSORS #####################################

             dataNames    = function() { colnames(_data)                       },
             targets      = function() { _targetVars                           },
             countMissing = function() { _missCounts                           },
             
###---------------------------------------------------------------------------###
             
             getControl = function () {
                 "Return the 'control' list"
                 list(convThresh        <- _convThresh
                      usePcStarts       <- _usePcStarts
                      center            <- _center
                      scale             <- _scale
                      adaptScales       <- _adaptScales
                      simpleIntercept   <- _simpleIntercept 
                      minPredCor        <- _minPredCor
                      miceIters         <- _miceIters
                      miceRidge         <- _miceRidge
                      miceMethod        <- _miceMethod
                      fimlStarts        <- _fimlStarts
                      preserveStructure <- _preserveStructure
                      optTraceLevel     <- _optTraceLevel
                      optCheckKkt       <- _optCheckKkt
                      optMethod         <- _optMethod
                      optBoundLambda    <- _optBoundLambda)
             },
             
###---------------------------------------------------------------------------###
             
             getImpDataset = function() {
                 "Fill missing values to produce a single imputed dataset"
                 data <- _data # Make a local copy of '_data'
                 
                 ## Randomly choose a posterior draw to use as imputations:
                 impRow       <-  sample(_impRowsPool, 1)
                 _impRowsPool <<- setdiff(_impRowsPool, impRow)
                 
                 for(j in _targetVars) {
                     impSam <- _gibbsOut[[j]]$imps[impRow, ]
                     data[_missList[[j]], j] <- impSam + _dataMeans[j]
                 }
                 
                 ## Restructure imputed data to match the raw data layout:
                 if(_preserveStructure)
                     data.frame(data, _ignoredColumns)[ , _rawNames]
             },

############################# DESCRIPTIVE FUNCTIONS #############################

             nVar     = function() { ncol(_data)                               },
             nObs     = function() { nrow(_data)                               },
             nPreds   = function() { nVar() - 1                                },
             nTargets = function() { length(_targetVars)                       },
             
########################## COMPLEX METHODS/SUBROUTINES ##########################
             
             applyMissCode = function() {
                 "Construct an integer-valued missing data code"
                 if(is.na(_missCode)) {
                     if(max(abs(_data), na.rm = TRUE) < 1.0) {
                         _missCode <<- -9
                     }
                     else {
                         codeMag   <-
                             floor(log10(max(abs(_data), na.rm = TRUE))) + 2
                         _missCode <<- -(10^codeMag - 1)
                     }
                 }
                 _data[is.na(_data)] <- _missCode
             },

###---------------------------------------------------------------------------###
             
             checkInputs = function() {
                 "Check the user inputs and isolate a set of target variables"
                 
                 ## Check for target variables. When no targets are given, all
                 ## incomplete variables not listed in 'ignoreVars' are imputed.
                 if(is.null(_targetVars)) {
                     if(_doImp) {
                         targetCandidates <-
                             colnames(_data)[!colnames(_data) %in% _ignoreVars]
                         warning("You did not specify any target variables, so I will impute the missing data on\nevery variable in 'data' that is not listed in 'ignoreVars'.\n")        
                     } else {
                         stop("Please specify a DV.")
                     }
                 } else {
                     targetCandidates <- _targetVars
                 }
                 
                 ## Make sure '_data' contains missing data that we can find:
                 if(is.na(_missCode)) {
                     rMat <- is.na(_data)
                 } else {
                     rMat <- _data == _missCode
                     
                     if(!any(rMat, na.rm = TRUE))
                         stop(paste0("The value you provided for 'missCode' (i.e., ",
                                     _missCode,
                                     ") does not appear anywhere in 'data'.\nAre you sure that ",
                                     _missCode,
                                     " encodes your missing data?\n")
                              )
                 }
                 
                 if(length(targetCandidates) > 1) 
                     completeTargets <- colMeans(rMat[ , targetCandidates]) == 0
                 else 
                     completeTargets <- mean(rMat[ , targetCandidates]) == 0
                 
                 if(_doImp & all(completeTargets)) 
                     stop("Your target variables appear to be fully observed. Did you forget to provide a\nvalue for 'missCode'?\n")
                 
                 ## Select the final set of target variables:
                 if(_doImp) {
                     _targetVars <<- targetCandidates[!completeTargets]
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
             
             scaleData = function(revert = FALSE) {
                 "Standardize the columns of _data"
                 nObs <- nrow(_data)
                 nVar <- ncol(_data)
                 
                 if(!revert) {# Doing initial scaling
                     if(_scale)
                         _dataScales <<- unlist(lapply(_data, sd))
                     else
                         _dataScales <<- rep(1, nVar)
                     
                     ## Mean center data:
                     if(_center) {
                         _dataMeans <<- colMeans(_data)
                         _data      <<- as.data.frame(
                             scale(_data, center = TRUE, scale = FALSE)
                         )
                     } else {
                         _dataMeans <<- rep(0, nVar)
                     }
                     
                     names(_dataMeans) <<- names(_dataScales) <- colnames(_data)
                     
                 } else {# Reverting the data to its original scaling
                     _data <<-
                         _data + data.frame(
                                     matrix(_dataMeans, nObs, nVar, byrow = TRUE)
                                 )
                 }
             },

###---------------------------------------------------------------------------###

             nameOutput = function() {
                 "Give the Gibb's sampling output pretty names"
                 if(_returnConvInfo) names(_rHatList) <- _targetVars
                 
                 if(_returnParams)
                     for(v in _targetVars) {
                         tmp <- setdiff(colnames(_data), v)
                         
                         colnames(_gibbsOut[[v]]$beta) <<- c("intercept", tmp)
                         colnames(_gibbsOut[[v]]$tau)  <<- tmp
                     }
             },


###---------------------------------------------------------------------------###

             calcRHat = function(sims) {
                 "Compute a single split-chain Potential Scale Reduction Factor"
                 subChainLen <- floor(length(sims) / 2)
                 nSubChains  <- _nChains * 2
                 
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
                 for(j in _targetVars) {
                     _rHats[[j]]$beta  <- apply(_gibbsOut[[j]]$beta, 2, calcRHat)
                     _rHats[[j]]$tau   <- apply(_gibbsOut[[j]]$tau,  2, calcRHat)
                     _rHats[[j]]$sigma <- calcRHat(gibbsOut[[j]]$sigma)
                 }
             }

###---------------------------------------------------------------------------###
             
             checkGibbsConv = function() {
                 "Check that the Gibb's sampler has converged everywhere"
                 for(j in _targetVars) {
                     ## Find nonconvergent Gibbs samples:
                     badBetaCount <- sum(_rHats[[j]]$beta > _convThresh)
                     maxBetaRHat  <- max(_rHats[[j]]$beta)
                     badTauCount  <- sum(_rHats[[j]]$tau > _convThresh)
                     maxTauRHat   <- max(__rHats[[j]]$tau)
                     badSigmaFlag <- _rHats[[j]]$sigma > _convThresh
                     
                     ## Return warnings about possible failures of convergence:
                     if(badBetaCount > 0) {
                         warning(paste0("While imputing ",
                                        j,
                                        ", Beta's final Gibbs sample may not have converged.\n",
                                        badBetaCount,
                                        " R-Hats > ",
                                        _convThresh,
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
                                        _convThresh,
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
                 }# CLOSE for(j in _targetVars)
             },
             
###---------------------------------------------------------------------------###
             
             startParams = function(restart = FALSE) {    
                 "Provide starting values for all parameters"

                 nTargets <- length(_targetVars)
                 if(restart) {
                     for(j in 1 : nTargets) {
                         _sigmaStarts[j]   <- mean(_gibbsOut[[j]]$sigma)
                         _tauStarts[ , j]  <- colMeans(_gibbsOut[[j]]$tau)
                         _betaStarts[ , j] <-
                             colMeans(_gibbsOut[[j]]$beta[ , -1])
                     }
                     return
                 }
                 
                 nRows    <- nrow(_data)
                 nObsVec  <- colSums(!is.na(_data))
                 nPreds   <- ncol(_data) - 1
                 
                 ## NOTE: We don't need to start the intercept. It's initial
                 ##       value will be sampled in the first iteration of the
                 ##       Gibbs sampler.
                 
                 ## Populate the starting values for Lambda:
                 if(!doBl) {
                     _lambdaMat <<- cbind(
                         matrix(_lambda1Starts, nTargets, 1),
                         matrix(_lambda2Starts, nTargets, 1)
                     )
                 } else {
                     _lambdaMat <<- cbind(matrix(_lambda1Starts, nTargets, 1), 0)
                 }# END if(!doBl)
                 
                 ## Populate starting values for betas, taus, and sigma:
                 _sigmaStarts <<- _dataScales[_targetVars]
                 
                 for(j in 1 : nTargets) {
                     if(!_doBl) {
                         lam1 <- _lambdaMat[j, 1]
                         lam2 <- _lambdaMat[j, 2]
                         
                         tauPriorScale <- (8 * lam2 * _sigmaStarts[j]) / lam1^2
                         
                         for(k in 1 : nPreds) {
                             tauDraw <- 0.0
                             while(tauDraw < 1.0)
                                 tauDraw <- rgamma(n     = 1,
                                                   shape = 0.5,
                                                   scale = tauPriorScale)
                             _tauStarts[k, j] <<- tauDraw
                         }
                         
                         betaPriorCov <- diag(
                             1 / ((lam2 / _sigmaStarts[j]) *
                                  (_tauStarts[ , j] / (_tauStarts[ , j] - 1.0))
                             )
                         )
                         
                         _betaStarts[ , j] <<-
                             rmvnorm(1, rep(0, nPreds), betaPriorCov)
                     } else {# We're doing BL
                         lam <- _lambdaMat[j, 1]
                         
                         _tauStarts[ , j] <<- rexp(nPreds, rate = (0.5 * lam^2))
                         
                         betaPriorCov <- _sigmaStarts[j] * diag(_tauStarts[ , j])
                         _betaStarts[ , j] <<-
                             rmvnorm(1, rep(0, nPreds), betaPriorCov)
                     }
                 }# CLOSE for(j in 1 : nTargets)
             },

###---------------------------------------------------------------------------###
             
             smoothLambda = function() {
                 i     <- iterations[1]
                 range <- (i - _smoothingWindow + 1) : i

                 for(j in _targetVars)
                     lambdaMat[j, ] <- colMeans(_lambdaHistory[[j]][range, ])
             }
             
         )# END MibrrFit$methods()
