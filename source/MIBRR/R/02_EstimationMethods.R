### Title:    Optimization and Gibbs Sampling Methods for MIBRR
### Author:   Kyle M. Lang
### Created:  2017-SEP-30
### Modified: 2019-FEB-01
### Notes:    This file will add optimization and Gibbs sampling methods to the
###           MibrrFit class.

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


MibrrFit$methods(

             doGibbs = function(phase = 1) {
                 "Run the Gibbs sampler to update the imputation model parameters"

                 respCounts <- nObs - missCounts

                 ## Get a new vector of seeds for the C++ samplers:
                 seedVec <- rep(NA, nTargets)
                 for(v in 1 : nTargets) {
                     sName <- paste0("mibrrStream", v)
                     .lec.ResetNextSubstream(sName)
                     seedVec[v] <- .lec.GetState(sName)[1]
                 }
                
                 gibbsOut <<-
                     runGibbs(data        = as.matrix(data),
                              nTargets    = nTargets,
                              missList    = missList[targetVars],
                              respCounts  = respCounts[targetVars],
                              lambda1     = lambdaMat[ , 1], 
                              lambda2     = lambdaMat[ , 2], # Ignored for BL
                              l1Parms     = l1Pars, # Ignored when
                              l2Parms     = l2Pars, # doMcem = TRUE
                              sigmaStarts = sigmaStarts,
                              tauStarts   = tauStarts,
                              betaStarts  = betaStarts,
                              burnSams    = sampleSizes[[phase]][1],
                              totalSams   = sum(sampleSizes[[phase]]),
                              penType     = penalty,
                              ridge       = ridge,
                              verbose     = verbose,
                              fullBayes   = !doMcem,
                              noMiss      = all(missCounts == 0),
                              savePpSams  = savePpSams,
                              seeds       = seedVec)
                 
                 names(gibbsOut) <<- targetVars
                 
                 ## Update the parameters' starting values:
                 if(doMcem) startParams(restart = TRUE)
             },
             
             eNetLL = function(lambdaVec, index) {
                 "Conditional loglikelihood function of Lambda (used for MCEM)"
                 
                 l1 <- lambdaVec[1]
                 l2 <- lambdaVec[2]
    
                 taus   <- gibbsOut[[index]]$tau
                 sigmas <- gibbsOut[[index]]$sigma
                 betas  <- gibbsOut[[index]]$beta
                 
                 p <- ncol(taus)
                 
                 e1 <- mean(
                     log(pgamma(l1^2 / (8 * sigmas * l2),
                                0.5,
                                lower.tail = FALSE) *
                         gamma(0.5)
                         )
                 )
                 e2 <- mean(
                     rowSums((taus / (taus - 1)) * betas[ , -1]^2) / sigmas
                 )
                 e3 <- mean(rowSums(taus) / sigmas)

                 ## Return the LL:
                 p * log(l1) - p * e1 - (l2 / 2) * e2 - (l1^2 / (8 * l2)) * e3
             },
             
             eNetGrad = function(lambdaVec, index) {
                 "Gradient function for the conditional LL of Lambda"
                 l1 <- lambdaVec[1]
                 l2 <- lambdaVec[2]
                 
                 taus   <- gibbsOut[[index]]$tau
                 sigmas <- gibbsOut[[index]]$sigma
                 betas  <- gibbsOut[[index]]$beta
                 
                 p   <- ncol(taus)
                 tmp <- l1^2 / (8 * sigmas * l2)
                 
                 e1 <- mean(
                 (1 / (pgamma(tmp, 0.5, lower.tail = FALSE) * gamma(0.5))) *
                 (1 / (sqrt(tmp) * exp(tmp))) * (1 / sigmas)
                 )
                 e2 <- mean(
                     rowSums((taus / (taus - 1)) * betas[ , -1]^2) / sigmas
                 )
                 e3 <- mean(rowSums(taus) / sigmas)
                 
                 w1 <- l1 / (4 * l2)
                 w2 <- l1^2 / (8 * l2^2)
                 
                 c((p / l1) + (p * w1 * e1) - (w1 * e3),  # dLL / dl1
                 (-p * w2 * e1) - (0.5 * e2) + (w2 * e3)) # dLL / dl2
             },             
                      
             optWrap = function(index, method, lowBounds) {
                 "Wrapper to allow optimx to run within lapply()"
                 optOut <- optimx(par     = lambdaMat[index, ],
                                  fn      = .self$eNetLL,
                                  gr      = .self$eNetGrad,
                                  method  = method,
                                  lower   = lowBounds,
                                  control = list(
                                      trace     = optTraceLevel,
                                      maximize  = TRUE,
                                      kkt       = optCheckKkt,
                                      follow.on = length(method) > 1
                                  ),
                                  index   = index)
                 
                 if(length(method) > 1) optOut <- optOut[nrow(optOut), ]
                 
                 ## Store the optimized lambdas:
                 lambdaMat[index, ] <<- coef(optOut)
                 
                 ## Return convergence code and KKT optimality checks
                 optOut[c("convcode", "kkt1", "kkt2")]
             },
             
             updateBlLambda = function(index) {
                 "Optimize lambda for the BL using the rule given in Park & Casella (2008)"
                 taus <- gibbsOut[[index]]$tau
                 p    <- ncol(taus)
                 
                 lambdaMat[index, 1] <<- sqrt((2 * p) / sum(colMeans(taus)))
             },

             optimizeLambda = function(iter) {
                 "Optimize the BEN or BL penalty parameters"
                 
                 ## Use simple update rule and return early when doing BL:
                 if(penalty == 1) {
                     lapply(1 : nTargets, .self$updateBlLambda)
                 }
                 else {
                     if(optBoundLambda) lowBounds <- c(1e-5, 1e-5)
                     else               lowBounds <- -Inf
                     
                     options(warn = ifelse(verbose, 0, -1))
                     
                     ## Define a location in which to sink unwanted output:
                     if(.Platform$OS.type == "unix") nullFile <- "/dev/null"
                     else                            nullFile <- "nul"
                     
                     if(optTraceLevel == 0) sink(nullFile) # Suppress optimx output
                     
                     ## Apply over targets to optimize lambdas:
                     convList <- lapply(X         = 1 : nTargets,
                                        FUN       = .self$optWrap,
                                        method    = optMethod,
                                        lowBounds = lowBounds)
                     
                     if(optTraceLevel == 0) sink()
                     options(warn = 0)
                 }# CLOSE if(penalty == 1); else
                 
                 for(j in 1 : nTargets) {
                     lambdaHistory[[j]][iter, ] <<- lambdaMat[j, ]
                     lambdaConv[[j]][iter, ]    <<- convList[[j]]
                     
                     ## Smooth Lambda estimates if beginning 'tuning' phase:
                     if(iter == iterations[1] & smoothingWindow > 1) {
                         smoothRange    <- (iter - smoothingWindow + 1) : iter         
                         lambdaMat[j, ] <<-
                             colMeans(lambdaHistory[[j]][smoothRange, ])        
                     }
                 }
             }
             
         )# END MibrrFit$methods()
