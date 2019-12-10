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
                              chain  = "integer",
                              data              = "data.frame",
                              targetVars        = "character",
                                        #ignoreVars        = "character",
                              iterations        = "integer",
                              sampleSizes       = "list",
                                        #missCode          = "integer",
                              seed              = "ANY",
                              doImp             = "logical",
                              doMcem            = "logical",
                                        #checkConv         = "logical",
                                        #verbose           = "logical",
                                        #convThresh        = "numeric",
                                        #lambda1Starts     = "numeric",
                                        #lambda2Starts     = "numeric",
                              l1Pars            = "numeric",
                              l2Pars            = "numeric",
                                        #usePcStarts       = "logical",
                              smoothingWindow   = "integer",
                                        #minPredCor        = "numeric",
                                        #miceIters         = "integer",
                                        #miceRidge         = "numeric",
                                        #miceMethod        = "character",
                                        #preserveStructure = "logical",
                              optTraceLevel     = "integer",
                              optCheckKkt       = "logical",
                              optMethod         = "character",
                              optBoundLambda    = "logical",
                                        #gibbsOut          = "list",
                                        #ignoredColumns    = "data.frame",
                                        #rawNames          = "character",
                                        #impRowsPool       = "integer",
                              missList          = "list",
                                        #nChains           = "integer",
                                        #rHats             = "list",
                                        #lambdaMat         = "matrix",
                                        #lambdaHistory     = "list",
                                        #lambdaConv        = "list",
                                        #betaStarts        = "matrix",
                                        #tauStarts         = "matrix",
                                        #sigmaStarts       = "numeric",
                                        #userMissCode      = "logical",
                                        #missCounts        = "integer",
                                        #nTargets          = "integer",
                                        #nVar              = "integer",
                                        #nPreds            = "integer",
                                        #nObs              = "integer",
                              totalIters        = "integer",
                              rng0              = "character",
                              userRng           = "character",
                              ridge             = "numeric",
                              penalty           = "integer",
                              savePpSams        = "logical",
                              useBetaMeans      = "logical",
                              optMaxRestarts    = "integer",
                              optRestartRatio   = "numeric",
                              optStrict         = "logical",
                              centerType        = "character",
                              dumpParamHistory  = "logical",
                              phHistoryLength   = "integer",
                              parameters        = "list"
                          )
                          )

###--------------------------------------------------------------------------###

MibrrChain$methods(
               
################################ CONSTRUCTOR ###################################
               
               initialize =
                   function(data        = data.frame(NULL),
                            targetVars  = "",
                                        #ignoreVars  = "",
                            iterations  = as.integer(NA),
                            sampleSizes = list(),
                            missCode    = as.integer(NA),
                                        #doImp       = as.logical(NA),
                            doMcem      = as.logical(NA),
                            verbose     = as.logical(NA),
                            seed        = NULL,
                            userRng     = "",
                            ridge       = 0.0,
                            penalty     = 2L,
                                        #nChains     = 1L
                            )
                   {
                       "Initialize an object of class MibrrFit"
                       data              <<- data
                       targetVars        <<- targetVars
                                        #ignoreVars        <<- ignoreVars
                       iterations        <<- iterations
                       sampleSizes       <<- sampleSizes
                       missCode          <<- missCode
                                        #doImp             <<- doImp
                       doMcem            <<- doMcem
                                        #checkConv         <<- TRUE
                       verbose           <<- verbose
                                        #convThresh        <<- 1.1
                                        #usePcStarts       <<- FALSE
                                        #minPredCor        <<- 0.3
                                        #miceIters         <<- 10L
                                        #miceRidge         <<- 1e-4
                                        #miceMethod        <<- "pmm"
                                        #preserveStructure <<- TRUE
                       optTraceLevel     <<- 0L
                       optCheckKkt       <<- TRUE
                       optMethod         <<- "L-BFGS-B"
                       optBoundLambda    <<- TRUE
                                        #nChains           <<- nChains
                       seed              <<- seed
                       userRng           <<- userRng
                       ridge             <<- ridge
                       penalty           <<- penalty
                       savePpSams        <<- FALSE
                       useBetaMeans      <<- FALSE
                       optMaxRestarts    <<- 5L
                       optRestartRatio   <<- 0.1
                       optStrict         <<- TRUE
                       centerType        <<- "median"
                       dumpParamHistory  <<- FALSE
                       phHistoryLength   <<- 10L

                       ## Initialize the parameter samples:
                       for(j in targetVars) 
                           parameters[[j]] <<-
                               MibrrSamples(target      = j,
                                            iterations  = sum(iterations),
                                            doMcem      = doMcem,
                                            penalty     = penalty,
                                            predVars    =
                                                setdiff(colnames(data), j),
                                            targetScale =
                                                sd(data[[j]], na.rm = TRUE)
                                            )
                   },
               
###--------------------------------------------------------------------------###
               
               prepStarts = function() {
                   "Prepare the Gibbs sampler's starting values"
                   
                   lambda1 <- lambda2 <- rep(NA, nTargets)
                   beta    <- matrix(NA, nPreds + 1, nTargets)
                   tau     <- matrix(NA, nPreds, nTargets)
                   sigma   <- rep(NA, nTargets)
                   
                   for(j in targetVars)
                       with(parameters[[j]],
                       {
                           lambda1[j] <- getLambda1()
                           lambda2[j] <- getLambda2()
                           beta[ , j] <- starts$beta
                           tau[ , j]  <- starts$tau
                           sigma[j]   <- starts$sigma
                       }
                       )
               },
               
###--------------------------------------------------------------------------###
               
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

                   ## Extract starting values from the 'parameters' field:
                   starts <- prepStarts()
                   
                   gibbsOut <<-
                       runGibbs(data         = as.matrix(data),
                                nTargets     = nTargets,
                                missList     = missList[targetVars],
                                respCounts   = respCounts[targetVars],
                                lambda1      = starts$lambda1,               #lambdaMat[ , 1], 
                                lambda2      = starts$lambda2,               #lambdaMat[ , 2], # Ignored for BL
                                l1Parms      = l1Pars, # Ignored when
                                l2Parms      = l2Pars, # doMcem = TRUE
                                sigmaStarts  = starts$sigma,                  #sigmaStarts,
                                tauStarts    = starts$tau,                    #tauStarts,
                                betaStarts   = starts$beta,                   #betaStarts,
                                burnSams     = sampleSizes[[phase]][1],
                                totalSams    = sum(sampleSizes[[phase]]),
                                penType      = penalty,
                                ridge        = ridge,
                                verbose      = verbose,
                                fullBayes    = !doMcem,
                                noMiss       = all(missCounts == 0),
                                savePpSams   = savePpSams,
                                useBetaMeans = useBetaMeans,
                                finalRep     = phase == 3,
                                seeds        = seedVec)
                   
                   names(gibbsOut) <<- targetVars

                   for(j in targetVars) {
                       ## Update the 'parameters' field:
                       parameters$setSamples(gibbsOut)
                       
                       ## Update the parameters' starting values:
                       if(doMcem)
                           parameters[[j]]$startParams(restart = TRUE)
                   }
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
                               ## Use the first set of estimates when proceeding:
                               parameters[[target]]$setLambdas(lamHat0)
                               tryOpt <-  FALSE
                               warning(msg, call. = FALSE, immediate. = TRUE)
                           }
                       }
                   }# CLOSE while(tryOpt)
                   
                   ## Return convergence info:
                   conv
               },

###--------------------------------------------------------------------------###

               updateBlLambda = function(target) {
                   "Optimize lambda for the BL using the rule given in Park & Casella (2008)"
                   taus <- parameters[[target]]$tau
                   p    <- ncol(taus)
                   
                   parameters[[target]]$setLambda1(
                                            sqrt((2 * p) / sum(colMeans(taus)))
                                        )
               },

###--------------------------------------------------------------------------###

               optimizeLambda = function(iter) {
                   "Optimize the BEN or BL penalty parameters"
                   
                   ## Use simple update rule and return early when doing BL:
                   if(penalty == 1) {
                       lapply(targetVars, .self$updateBlLambda)
                   }
                   else {
                       if(optBoundLambda) lowBounds <- c(1e-5, 1e-5)
                       else               lowBounds <- -Inf
                       
                                        #options(warn = ifelse(verbose, 0, -1))
                       
                       ## Define a location in which to sink unwanted output:
                       if(.Platform$OS.type == "unix") nullFile <- "/dev/null"
                       else                            nullFile <- "nul"
                       
                       if(optTraceLevel == 0) sink(nullFile) # Suppress optimx output
                       
                       ## Apply over targets to optimize lambdas:
                       convList <- lapply(X         = targetVars,
                                          FUN       = .self$optWrap,
                                          method    = optMethod,
                                          lowBounds = lowBounds)
                       
                       if(optTraceLevel == 0) sink()
                       options(warn = 0)
                   }# CLOSE if(penalty == 1); else
                   
                   for(j in 1 : nTargets) {
                                        #lambdaHistory[[j]][iter, ] <<- lambdaMat[j, ]

                       if(penalty == 2)
                           parameters[[target]]$setLambdaConv(convList[[j]])
                       
                       ## Smooth Lambda estimates if beginning 'tuning' phase:
                                        #if(iter == iterations[1] & smoothingWindow > 1) {
                                        #    smoothRange    <- (iter - smoothingWindow + 1) : iter         
                                        #    lambdaMat[j, ] <<-
                                        #        colMeans(lambdaHistory[[j]][smoothRange, ])        
                   }
               }
}

)# END MibrrFit$methods()
