### Title:    MibrrSamples Reference Class Definition
### Author:   Kyle M. Lang
### Created:  2017-12-09
### Modified: 2019-12-09
### Note:     The MibrrSamples class holds the parameter samples for one target
###           variable and one Markov chain

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

MibrrSamples <- setRefClass("MibrrSamples",
                            fields = list(
                                iter       = "integer",
                                target     = "character",
                                predVars   = "character",
                                nIters     = "integer",
                                targetSd   = "numeric",
                                beta       = "matrix",
                                tau        = "matrix",
                                sigma      = "numeric",
                                lambda1    = "numeric",
                                lambda2    = "numeric",
                                imps       = "matrix",
                                ppSams     = "matrix",
                                starts     = "list",
                                lambdaConv = "data.frame",
                                doMcem     = "logical",
                                penalty    = "integer",
                                centerType = "character"
                            )
                            )

###--------------------------------------------------------------------------###

MibrrSamples$methods(
                 
################################ CONSTRUCTOR ###################################
                 
                 initialize = function(target   = "",
                                       predVars = "",
                                       targetSd = as.numeric(NA),
                                       nIters   = as.integer(NA),
                                       doMcem   = as.logical(NA),
                                       penalty  = as.integer(NA)
                                       )
                 {
                     "Initialize an object of class MibrrSamples"
                     target   <<- target
                     predVars <<- predVars
                     targetSd <<- targetSd
                     nIters   <<- nIters
                     doMcem   <<- doMcem
                     penalty  <<- penalty
                     iter     <<- 1L
                     
                     ## Initialize penalty parameter-related stuff:
                     starts$lambda1 <<- 0.5
                     starts$lambda2 <<- length(predVars) / 10
                     
                     if(doMcem) {
                         lambda1 <<- vector("numeric", nIters)
                         lambda2 <<- vector("numeric", nIters)
                         
                         lambdaConv <<-
                             data.frame(code = vector("integer", nIters),
                                        kkt1 = vector("logical", nIters),
                                        kkt2 = vector("logical", nIters)
                                        )
                     }
                     
                     ## Provide starting values for the parameters:
                     .self$startParams()
                 },
################################## ACCESSORS ###################################

                 getLambda1 = function() lambda1[iter],
                 
###--------------------------------------------------------------------------###

                 getLambda2 = function() lambda2[iter],

###--------------------------------------------------------------------------###

                 getLambdas = function() c(lambda1[iter], lambda2[iter]),
                 
################################### MUTATORS ###################################
                 
                                        #setControl = function(x) {
                                        #    "Assign the control parameters"
                                        #    
                                        #    ## Get the fields for each class:
                                        #    fields <- getRefClass(class(.self))$fields()
                                        #    
                                        #    ## Assign the control list entries to the correct classes:
                                        #    for(n in names(x))
                                        #        if(n %in% names(fields))
                                        #            field(n, cast(x[n], fields[n]))
                                        #},
                 
###--------------------------------------------------------------------------###
                 
                 incIter = function() iter <<- iter + 1L,
                 
###--------------------------------------------------------------------------###
                 
                 setSamples = function(gibbsOut) {
                     ## Update the primary parameter samples:
                     for(x in c("beta", "tau", "sigma"))
                         field(x, gibbsOut[[target]][[x]])
                     
                     ## Update the imputations:
                     if(length(gibbsOut[[target]]$imps) > 1)
                         field("imps", gibbsOut[[target]]$imps)
                     
                     ## Update the posterior predictive samples:
                     if(length(gibbsOut[[target]]$ppSams) > 1)
                         field("ppSams", gibbsOut[[target]]$ppSams)
                     
                     ## Update the penalty parameters' samples:
                     if(!doMcem & penalty != 0) {
                         field("lambda1", gibbsOut[[target]]$lambda[ , 1])
                         
                         if(penalty == 2)
                             field("lambda2", gibbsOut[[target]]$lambda[ , 2])
                     }
                 },
                 
###--------------------------------------------------------------------------###
                 
                 setLambda1 = function(l1) lambda1[iter] <<- l1,
                 
###--------------------------------------------------------------------------###

                 setLambdas = function(lams) {
                     lambda1[iter] <<- lams[1]
                     lambda2[iter] <<- lams[2]
                 },

###--------------------------------------------------------------------------###
                 
                 setLambdaConv = function(x) lambdaConv[iter, ] <<- x,
                 
###--------------------------------------------------------------------------###
                 
                 startParams = function(restart = FALSE) {    
                     "Provide starting values for all parameters"
                     
                     if(restart) {
                         ## Choose the type of central tendency to use:
                         cenTen <- switch(centerType,
                                          mean   = mean,
                                          median = median,
                                          mode   = numMode)
                         
                         starts$sigma <<- cenTen(sigma)
                         starts$tau   <<- apply(tau, 2, cenTen)
                         starts$beta  <<- apply(beta, 2, cenTen)

                         starts$lambda1 <<- getLambda1()
                         starts$lambda2 <<- getLambda2()
                         
                         return()
                     }
                     
                     ## NOTE: We don't need real starting values for the
                     ##       intercepts. Their initial values will be sampled
                     ##       in the first iteration of the Gibbs sampler.
                     
                     ## Populate the starting values for Lambda:
                     lambda1[1] <<- starts$lambda1
                     lambda2[1] <<- starts$lambda2
                     
                     ## Populate starting values for betas, taus, and sigma:
                     shape        <-  targetSd^2 / (0.1 * targetSd)
                     scale        <-  (0.1 * targetSd) / targetSd
                     starts$sigma <<- rgamma(1, shape = shape, scale = scale)

                     nPreds <- length(predVars)
                     
                     if(penalty == 2) {     # We're doing BEN
                         lam1 <- lambda1[1]
                         lam2 <- lambda2[1]
                         
                         tauPriorScale <- (8 * lam2 * starts$sigma) / lam1^2

                         starts$tau <<- rep(NA, nPreds)
                         for(k in 1 : nPreds) {
                             tauDraw <- 0.0
                             while(tauDraw < 1.0)
                                 tauDraw <- rgamma(n     = 1,
                                                   shape = 0.5,
                                                   scale = tauPriorScale)
                             starts$tau[k] <<- tauDraw
                         }
                         
                         betaPriorCov <-
                             with(starts,
                                  diag(
                                      1 / ((lam2 / sigma) * (tau / (tau - 1.0)))
                                  )
                                  )
                     }
                     else if(penalty == 1) {# We're doing BL
                         lam          <-  lambda1[1]
                         starts$tau   <<- rexp(nPreds, rate = (0.5 * lam^2))
                         betaPriorCov <-  with(starts, sigma * diag(tau))
                     }
                     else                   # We're doing basic ridge 
                         betaPriorCov <- diag(rep(starts$sigma, nPreds)) 
                     
                     starts$beta <<-
                         c(0, rmvnorm(1, rep(0, nPreds), betaPriorCov))
                 }
                 
             )# END MibrrSamples$methods()

