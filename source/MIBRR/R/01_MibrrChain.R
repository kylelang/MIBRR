### Title:    MibrrChain Reference Class Definition
### Author:   Kyle M. Lang
### Created:  2017-09-30
### Modified: 2019-12-10
### Notes:    The MibrrChain class hold the methods and metadata associated
###           with one Markov chain

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


MibrrChain <- setRefClass("MibrrChain",
                          fields = list(
                              chain           = "integer",
                              data            = "data.frame",
                              targetVars      = "character",
                              sampleSizes     = "list",
                              doMcem          = "logical",
                              verbose         = "logical",
                              l1Pars          = "numeric",
                              l2Pars          = "numeric",
                              missList        = "list",
                              nTargets        = "integer",
                              ridge           = "numeric",
                              penalty         = "integer",
                              savePpSams      = "logical",
                              useBetaMeans    = "logical",
                              optTraceLevel   = "integer",
                              optCheckKkt     = "logical",
                              optMethod       = "character",
                              optBoundLambda  = "logical",
                              optMaxRestarts  = "integer",
                              optRestartRatio = "numeric",
                              optStrict       = "logical",
                              parameters      = "list"
                          )
                          )

###--------------------------------------------------------------------------###

data.frame(NULL)

MibrrChain$methods(
               
################################ CONSTRUCTOR ###################################
               
               initialize =
                   function(chain       = as.integer(NA),
                            data        = data.frame(),
                            targetVars  = "",
                            iterations  = as.integer(NA),
                            sampleSizes = list(),
                            missList    = list(),
                            doMcem      = as.logical(NA),
                            verbose     = as.logical(NA),
                            ridge       = 0.0,
                            penalty     = 2L,
                            control     = list()
                            )
                   {
                       "Initialize an object of class MibrrFit"
                       chain       <<- chain
                       data        <<- data
                       targetVars  <<- targetVars
                       nTargets    <<- length(targetVars)
                       sampleSizes <<- sampleSizes
                       missList    <<- missList
                       doMcem      <<- doMcem
                       verbose     <<- verbose
                       ridge       <<- ridge
                       penalty     <<- penalty
                       
                       ## Initialize the parameter samples:
                       for(j in targetVars) {
                           parameters[[j]] <<-
                               MibrrSamples(target   = j,
                                            nIters   = iterations,
                                            doMcem   = doMcem,
                                            penalty  = penalty,
                                            targetSd = sd(data[[j]]),
                                            predVars =
                                                setdiff(colnames(data), j)
                                            )

                           ## Set control parameters for the 'MibrrSamples'
                           ## objects:
                           setControl(x = control, where = parameters[[j]])
                       }
                   },
                          
########################### ESTIMATION ROUTINES ################################
               
               prepStarts = function() {
                   "Prepare the Gibbs sampler's starting values"
                   
                   lambda1 <- lambda2 <- rep(NA, nTargets)
                   beta    <- matrix(NA, ncol(data), nTargets)
                   tau     <- matrix(NA, ncol(data) - 1, nTargets)
                   sigma   <- rep(NA, nTargets)
                   
                   colnames(beta) <- colnames(tau) <- names(sigma) <-
                       names(lambda1) <- names(lambda2) <- targetVars
                   
                   for(j in targetVars) {
                       lambda1[j] <- parameters[[j]]$starts$lambda1
                       lambda2[j] <- parameters[[j]]$starts$lambda2
                       beta[ , j] <- parameters[[j]]$starts$beta
                       tau[ , j]  <- parameters[[j]]$starts$tau
                       sigma[j]   <- parameters[[j]]$starts$sigma
                   }
                   
                   list(lambda1 = lambda1,
                        lambda2 = lambda2,
                        beta    = beta,
                        tau     = tau,
                        sigma   = sigma)
               },
               
###--------------------------------------------------------------------------###
               
               doGibbs = function(phase = 1) {
                   "Run the Gibbs sampler to update the imputation model parameters"
                   
                   respCounts <- nrow(data) - sapply(missList, length)

                   ## Get a new vector of seeds for the C++ samplers:
                   #seedVec        <- rep(NA, nTargets)
                   #names(seedVec) <- targetVars
                   #for(j in targetVars) {
                   #    sName <- paste0("c", chain, j)
                   #    .lec.ResetNextSubstream(sName)
                   #    seedVec[j] <- .lec.GetState(sName)[1]
                   #}

                   seedVec <- rep(235711, nTargets)
                   
                   ## Extract starting values from the 'parameters' field:
                   starts <- prepStarts()
                   
                   gibbsOut <-
                       MIBRR:::runGibbs(data         = as.matrix(data),
                                nTargets     = nTargets,
                                missList     = missList[targetVars],
                                respCounts   = respCounts[targetVars],
                                lambda1      = starts$lambda1, 
                                lambda2      = starts$lambda2, # Ignored for BL
                                l1Parms      = l1Pars, # Ignored when
                                l2Parms      = l2Pars, # doMcem = TRUE
                                sigmaStarts  = starts$sigma,
                                tauStarts    = starts$tau,
                                betaStarts   = starts$beta,
                                burnSams     = sampleSizes[[phase]][1],
                                totalSams    = sum(sampleSizes[[phase]]),
                                penType      = penalty,
                                ridge        = ridge,
                                verbose      = verbose,
                                fullBayes    = !doMcem,
                                noMiss       = all(respCounts == nrow(data)),
                                savePpSams   = savePpSams,
                                useBetaMeans = useBetaMeans,
                                finalRep     = phase == 3,
                                seeds        = seedVec)#,
                                #chain        = chain)
                   
                   names(gibbsOut) <- targetVars

                   ## Update the 'parameters' field:
                   for(j in targetVars)
                       parameters[[j]]$setSamples(gibbsOut)
               },
               
###--------------------------------------------------------------------------###

               eNetLL = function(lambdas, target) {
                   "Conditional loglikelihood function of Lambda (used for MCEM)"
                   
                   l1 <- lambdas[1]
                   l2 <- lambdas[2]
                   
                   taus   <- parameters[[target]]$tau
                   sigmas <- parameters[[target]]$sigma
                   betas  <- parameters[[target]]$beta
                   
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

###--------------------------------------------------------------------------###

               eNetGrad = function(lambdas, target) {
                   "Gradient function for the conditional LL of Lambda"
                   l1 <- lambdas[1]
                   l2 <- lambdas[2]
                   
                   taus   <- parameters[[target]]$tau
                   sigmas <- parameters[[target]]$sigma
                   betas  <- parameters[[target]]$beta
                   
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

###--------------------------------------------------------------------------###

               optWrap = function(target, method, lowBounds) {
                   "Wrapper to allow optimx to run within lapply()"

                   ## Store the current value of Lambda:
                   lam0 <- parameters[[target]]$getLambdas()
                   
                   rep    <- 0
                   tryOpt <- TRUE
                   while(tryOpt) {
                       ## Optimize Lambda:
                       optOut <- optimx(par     = lam0,
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
                                        target   = target)
                       
                       if(length(method) > 1) optOut <- optOut[nrow(optOut), ]
                       
                       ## Check convergence and KKT optimality:
                       conv  <- optOut[c("convcode", "kkt1", "kkt2")]
                       check <- conv[1] == 0 & conv[2] & conv[3]
                       
                       if(check) { 
                           ## Increment the MCEM iteration counter:
                           parameters[[target]]$incIter()
                           
                           ## Store the optimized lambdas:
                           parameters[[target]]$setLambdas(coef(optOut))
                           tryOpt <-  FALSE
                       }
                       else {
                           ## Freeze the first set of estimates:
                           if(rep == 0) lamHat0 <- coef(optOut)
                           
                           ## Add some noise to the starting values:
                           tmp  <- parameters[[target]]$getLambdas()
                           lam0 <- tmp + rnorm(2, 0, optRestartRatio * tmp)
                           rep  <- rep + 1
                           
                           ## Print a warning about the failure:
                           warning(paste0("I failed to optimize lambda for ",
                                          target,
                                          ". So, I'll try restart number ",
                                          rep,
                                          " of ",
                                          optMaxRestarts,
                                          ". New value = ", lam0, "."),
                                   call.      = FALSE,
                                   immediate. = TRUE
                                   )
                       }
                       
                       ## Signal a condition if we reach the maximum number of
                       ## restarts
                       if(rep == optMaxRestarts) {
                           msg <- paste0("Lambda for ",
                                         target,
                                         " could not be optimized after ",
                                         rep,
                                         " restarts.")
                           if(optStrict) {
                               sink() # Stop sinking output before throwing error
                               stop(msg, call. = FALSE)
                           }
                           else {
                               ## Increment the MCEM iteration counter:
                               parameters[[target]]$incIter()
                               
                               ## Use the first estimates when proceeding:
                               parameters[[target]]$setLambdas(lamHat0)
                               tryOpt <-  FALSE
                               warning(msg, call. = FALSE, immediate. = TRUE)
                           }
                       }
                   }# CLOSE while(tryOpt)
                   
                   ## Save convergence info:
                   parameters[[target]]$setLambdaConv(conv)
               },

###--------------------------------------------------------------------------###

               updateBlLambda = function(target) {
                   "Optimize lambda for the BL using the rule given in Park & Casella (2008)"
                   taus <- parameters[[target]]$tau
                   p    <- ncol(taus)
                   
                   ## Increment the MCEM iteration counter:
                   parameters[[target]]$incIter()

                   ## Save the updated lambda1 value:
                   parameters[[target]]$setLambda1(
                                            sqrt((2 * p) / sum(colMeans(taus)))
                                        )
               },

###--------------------------------------------------------------------------###

               optimizeLambda = function() {
                   "Optimize the BEN or BL penalty parameters"
                   
                   ## Use simple update rule and return early when doing BL:
                   if(penalty == 1) {
                       lapply(targetVars, .self$updateBlLambda)
                   }
                   else {
                       if(optBoundLambda) lowBounds <- c(1e-5, 1e-5)
                       else               lowBounds <- -Inf
                       
                       ## Define a location in which to sink unwanted output:
                       if(.Platform$OS.type == "unix") nullFile <- "/dev/null"
                       else                            nullFile <- "nul"
                       
                       if(optTraceLevel == 0) sink(nullFile) # Suppress optimx output
                       
                       ## Apply over targets to optimize lambdas:
                       lapply(X         = targetVars,
                              FUN       = .self$optWrap,
                              method    = optMethod,
                              lowBounds = lowBounds)
                       
                       if(optTraceLevel == 0) sink()
                   }# CLOSE if(penalty == 1); else

                   ## Update starting values for the next Gibbs sampler run:
                   for(j in targetVars)
                       parameters[[j]]$startParams(restart = TRUE)
                   
                   ## Smooth Lambda estimates if beginning 'tuning' phase:
                                        #if(iter == iterations[1] & smoothingWindow > 1) {
                                        #    smoothRange    <- (iter - smoothingWindow + 1) : iter         
                                        #    lambdaMat[j, ] <<-
                                        #        colMeans(lambdaHistory[[j]][smoothRange, ])        
                                        #}
               }
               
           )# END MibrrFit$methods()
