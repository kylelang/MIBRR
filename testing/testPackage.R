### Title:    Test MIBRR Package
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2017-OCT-25
### Purpose:  Script to help test the MIBRR package

rm(list = ls(all = TRUE))

library(mibrr)
library(mitools)
library(psych)

source("testingSupportFunctions.R")

data(bfi)
tmp <- na.omit(bfi)

ed.d           <- model.matrix(~factor(tmp$education))[ , -1]
colnames(ed.d) <-
    c("finish_hs", "some_college", "college_grad", "graduate_degree")

male            <- tmp$gender
male[male == 2] <- 0

cn   <- setdiff(colnames(bfi), c("gender", "education"))
bfi2 <- data.frame(tmp[ , cn], male, ed.d)

rownames(bfi2) <- NULL

targets  <- list(mcar = paste0("A", c(1 : 5)),
                 mar  = paste0("C", c(1 : 5)),
                 mnar = paste0("E", c(1 : 5))
                 )
pm       <- list(mcar = 0.2, mar = 0.3, mnar = 0.1)
snr      <- list(mar = 5, mnar = 5)
marPreds <- c("age",
              "male",
              "finish_hs",
              "some_college",
              "college_grad",
              "graduate_degree")

dat1 <- bfi2[sample(c(1 : nrow(bfi2)), 500), ]

dat2 <- imposeMissing(data    = bfi2,
                      targets = targets,
                      preds   = marPreds,
                      pm      = pm,
                      snr     = snr)$data 

mibenOut <- miben(data           = dat2,
                  nImps          = 25,
                  targetVars     = unlist(targets),
                  ignoreVars     = NULL,
                  iterations     = c(50, 10),
                  returnConvInfo = TRUE,
                  returnParams   = TRUE,
                  verbose        = TRUE)
ls(mibenOut)

mibenOut$lambdaHistory
mibenOut$rHats
mibenOut$imps[[1]]

miblOut <- mibl(data           = dat3,
                nImps          = 25,
                targetVars     = unlist(targets),
                ignoreVars     = NULL,
                iterations     = c(50, 10),
                returnConvInfo = TRUE,
                returnParams   = TRUE,
                verbose        = TRUE)

miceOut <- mice(dat2, m = 25, maxit = 5, method = "norm")



tmp2              <- complete(miceOut)
tmp2[is.na(tmp2)] <- -99999

debug(mibrr)
undebug(mibrr)

out1 <- mibrr(doBl           = FALSE,
              doImp          = TRUE,
              data           = tmp2,
              nImps          = 5,
              targetVars     = c("y", paste0("x", c(1 : 3))),
              ignoreVars     = "idNum",
              iterations     = c(20, 10),
              sampleSizes    = list(c(100, 50), c(1000, 500), c(10000, 5000)),
              missCode       = -99999,
              returnConvInfo = TRUE,
              returnParams   = TRUE,
              verbose        = TRUE,
              seed           = 235711,
              control        = list(optTraceLevel = 0)
              )

out1 <- miben(data           = tmp2,
              nImps          = 100,
              targetVars     = c("y", paste0("x", c(1 : 3))),
              ignoreVars     = "idNum",
              #iterations     = c(20, 10),
              #sampleSizes    = list(c(, 50), c(1000, 500), c(10000, 5000)),
              missCode       = -99999,
              returnConvInfo = TRUE,
              returnParams   = TRUE,
              verbose        = TRUE,
              seed           = 235711#,
              #control        = list(optTraceLevel = 0)
              )

ls(out1)

par(mfcol = c(2, 4))

for(j in out1$lambdaHistory) {
    plot(j[ , 1], type = "l")
    plot(j[ , 2], type = "l")
}

fits <- lapply(out1$imps, function(x) lm(y ~ x1 + x2 + x3, data = x))

tmp3 <- tmp2
tmp3[tmp3 == -99999] <- NA

out0 <- lm(y ~ x1 + x2 + x3, data = tmp3)
summary(out0)

MIcombine(fits)

ls(out1)

out1$lambdaHistory
out1$params
out1$imps

optOut <- readRDS("optOut.rds")

do.call(rbind, optOut)

control <- list(optCheckKkt = TRUE)


lapply(optOut, function(x) 

##### DEGUB OPTIMIZATION #####

### The conditional loglikelihood function of Lambda for use during the
### empirical bayes updating.
eNetLL <- function(lambdaVec, gibbsState) {
    l1 <- lambdaVec[1]
    l2 <- lambdaVec[2]
    
    taus   <- gibbsState$tau
    sigmas <- gibbsState$sigma
    betas  <- gibbsState$beta
    
    p <- ncol(taus)
    
    e1 <- mean(
        log(pgamma(l1^2 / (8 * sigmas * l2), 0.5, lower = FALSE) * gamma(0.5))
    )
    e2 <- mean(rowSums((taus / (taus - 1)) * betas[ , -1]^2) / sigmas)
    e3 <- mean(rowSums(taus) / sigmas)
    
    p * log(l1) - p * e1 - (l2 / 2) * e2 - (l1^2 / (8 * l2)) * e3 # LL
}# END eNetLL()


### The gradient function for the conditional LL of Lambda:
eNetGrad <- function(lambdaVec, gibbsState)
{
    l1 <- lambdaVec[1]
    l2 <- lambdaVec[2]

    taus   <- gibbsState$tau
    sigmas <- gibbsState$sigma
    betas  <- gibbsState$beta
    
    p   <- ncol(taus)
    tmp <- l1^2 / (8 * sigmas * l2)
    
    e1 <- mean(
    (1 / (pgamma(tmp, 0.5, lower = FALSE) * gamma(0.5))) *
    (1 / (sqrt(tmp) * exp(tmp))) * (1 / sigmas)
    )
    e2 <- mean(rowSums((taus / (taus - 1)) * betas[ , -1]^2) / sigmas)
    e3 <- mean(rowSums(taus) / sigmas)

    w1 <- l1 / (4 * l2)
    w2 <- l1^2 / (8 * l2^2)
    
    c((p / l1) + (p * w1 * e1) - (w1 * e3),  # dLL / dl1
    (-p * w2 * e1) - (0.5 * e2) + (w2 * e3)) # dLL / dl2
}# END eNetGrad()


## Wrapper to allow optimx to run within lapply():
optWrap <- function(targetIndex,
                    lambdaMat,
                    optFun,
                    optGrad,
                    optMethod,
                    optLower,
                    optHessian,
                    optControl,
                    myGibbs)
{
    optOut <- optimx(par        = lambdaMat[targetIndex, ],
                     fn         = optFun,
                     gr         = optGrad,
                     method     = optMethod,
                     lower      = optLower,
                     hessian    = optHessian,
                     control    = optControl,
                     gibbsState = myGibbs[[targetIndex]])
    
    if(length(optMethod) > 1) optOut <- optOut[nrow(optOut), ]
    
    tmpList <- list()
    if(optHessian) {
        hessMat      <- attr(optOut, "details")[ , "nhatend"][[1]]
        tmpList$vcov <- solve(-hessMat)
    }
    
    if(optControl$kkt) tmpList$kktFlags <- c(optOut$kkt1, optOut$kkt2)
    
    tmpList$lambda <- c(optOut[[1]], optOut[[2]])
    tmpList
}# END optWrap()


## Optimize the penalty parameters via numerical maximization of eNetLL():
optimizeLambda <- function(lambdaMat,
                           gibbsState,
                           printFlag  = TRUE,
                           returnVCOV = FALSE,
                           controlParms)
{
    optMethod <- controlParms$method
    useSeqOpt <- length(optMethod) > 1
    
    if(controlParms$boundLambda) {
        lowBounds <- c(0, 0)
        optMethod <- "L-BFGS-B"
    } else {
        lowBounds <- -Inf
    }
    
    options(warn = ifelse(controlParms$showWarns, 0, -1))
    
    if(!printFlag) sink("/dev/null")
    
    optList <- lapply(c(1 : nrow(lambdaMat)),
                      FUN        = optWrap,
                      lambdaMat  = lambdaMat,
                      optFun     = eNetLL,
                      optGrad    = eNetGrad,
                      optMethod  = optMethod,
                      optLower   = lowBounds,
                      optHessian = returnVCOV,
                      optControl =
                          list(trace     = controlParms$traceLevel,
                               maximize  = TRUE,
                               kkt       = controlParms$checkKKT,
                               follow.on = useSeqOpt),
                      myGibbs    = gibbsState)
    
    if(!printFlag) sink()
    options(warn = 0)
    
    optList
}# END optimizeLambda()


gibbsOut   <- readRDS("gibbsOut3.rds")
gibbsState <- gibbsOut[[1]]

grad(func = eNetLL, x = c(0.5, 1.5), gibbsState = gibbsState)
eNetGrad(c(0.5, 1.5), gibbsState)

ls(gibbsOut[[1]])


## Test MVN sampler:
nObs <- 500000
mvnMu <- rep(10, 3)
mvnSigma <- matrix(5, 3, 3)
diag(mvnSigma) <- 20

out1.1 <- mibrr::drawMVN(nObs, mvnMu, mvnSigma)
out1.2 <- rmvnorm(nObs, mvnMu, mvnSigma)

par(mfrow = c(1, 3))
for(i in 1 : length(mvnMu)) {
    plot(density(out1.1[ , i]), col = "red")
    lines(density(out1.2[ , i]), col = "blue")
}

## Test inverse gamma sampler:
gamShape <- 10
gamScale <- 10

out2.1 <- mibrr::drawInvGamma(nObs, gamShape, gamScale)
out2.2 <- rinvgamma(nObs, gamShape, gamScale)

plot(density(out2.1), col = "red")
lines(density(out2.2), col = "blue")

## Test inverse gaussian sampler:
igMu = 1
igLam = 2

out3.1 <- mibrr::drawInvGauss(nObs, igMu, igLam)
out3.2 <- rinvgauss(nObs, igMu, igLam)

plot(density(out3.1), col = "red")
lines(density(out3.2), col = "blue")

## Test the incomplte gamma calculation:
incGamShape <- 10
incGamCut <- 5

out4.1 <- mibrr::calcIncGamma(incGamShape, incGamCut, FALSE)
out4.2 <- pgamma(q = incGamCut,
                 shape = incGamShape,
                 lower = FALSE) * gamma(incGamShape)

out4.1 - out4.2

## Test MIBEN and MIBL:
data(mibrrExampleData)

miceOut <- mice(mibrrExampleData,
                m = 1,
                maxit = 100)

dat0 <- complete(miceOut)
dat1 <- data.frame(mibrrExampleData[ , 1 : 5], dat0[ , 6 : 17])

debug(miben)
undebug(miben)

testOut <- miben(rawData      = dat1,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                                        #iterations   = c(2, 1),
                                        #sampleSizes  = c(100, 500, 1000),
                 returnParams = TRUE,
                 verbose      = TRUE,
                 control      = list(simpleIntercept = TRUE,
                                     adaptScales     = TRUE)
                 )

fitOut <- lapply(testOut$imps,
                 FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                 )
MIcombine(fitOut)

summary(lm(y ~ x1 + x2 + x3, data = dat1))


testOut2 <- mibl(rawData      = mibrrExampleData,
                 targetVars   = c("y", paste0("x", c(1 : 3))),
                 ignoreVars   = "idNum",
                 returnParams = TRUE)

fitOut2 <- lapply(testOut2$imps,
                  FUN = function(x) lm(y ~ x1 + x2 + x3, data = x)
                  )
MIcombine(fitOut2)


## Test BEN and BL:
library(mibrr)

alpha  <- 0.5
nPreds <- 175

parms             <- list()
parms$nObs        <- 100
parms$nPreds      <- nPreds
parms$r2          <- 0.2
parms$collin      <- 0.3
parms$beta        <- matrix(c(alpha, runif(nPreds, 0.3, 0.6)))
parms$means       <- runif(nPreds, 0, 1)
parms$scales      <- rep(1, nPreds)
parms$center      <- TRUE
parms$scale       <- TRUE
parms$simpleInt   <- FALSE
parms$verbose     <- FALSE
parms$postN       <- 5000
parms$adaptScales <- TRUE

testFun(1, parms)

rp <- 1

colMeans(dat1)

testFun <- function(rp, parms) {
    nObs   <- parms$nObs
    nPreds <- parms$nPreds
    r2     <- parms$r2
    collin <- parms$collin
    beta   <- parms$beta
    means  <- parms$means
    scales <- parms$scales

    dat1 <- simulateData(nObs   = nObs,
                         nPreds = nPreds,
                         r2     = r2,
                         collin = collin,
                         beta   = beta,
                         means  = means,
                         scales = scales)
    
    testOut <- ben(rawData = dat1,
                   y       = "y",
                   X       = paste0("x", c(1 : nPreds)),
                   verbose = parms$verbose,
                   control =
                       list(center          = parms$center,
                            scale           = parms$scale,
                            adaptScales     = parms$adaptScales,
                            simpleIntercept = parms$simpleInt)
                   )
    
    testOut2 <- bl(rawData = dat1,
                   y       = "y",
                   X       = paste0("x", c(1 : nPreds)),
                   verbose = parms$verbose,
                   control =
                       list(center          = parms$center,
                            scale           = parms$scale,
                            adaptScales     = parms$adaptScales,
                            simpleIntercept = parms$simpleInt)
                   )
    
    testForm <-
        as.formula(paste0("y ~ ", paste0("x", c(1 : nPreds), collapse = " + ")))
    lmOut <- lm(testForm, data = dat1)
    
    nTests <- 100
    mseMat <- matrix(NA, nTests, 3)
    for(i in 1 : nTests) {
        testDat <- simulateData(nObs   = nObs,
                                nPreds = nPreds,
                                r2     = r2,
                                collin = collin,
                                beta   = beta,
                                means  = means,
                                scales = scales)
        
        predOut1 <- predictMibrr(object  = testOut,
                                 newData = as.matrix(testDat[ , -1])
                                 )
        
        predOut2 <- predictMibrr(object  = testOut2,
                                 newData = as.matrix(testDat[ , -1])
                                 )
        
        predOut3 <- predict(lmOut, newdata = testDat[ , -1])
        
        mseMat[i, 1] <- mean((predOut1 - dat1$y)^2)
        mseMat[i, 2] <- mean((predOut2 - dat1$y)^2)
        mseMat[i, 3] <- mean((predOut3 - dat1$y)^2)
    }
    
    outMat           <- rbind(colMeans(mseMat), apply(mseMat, 2, var))
    colnames(outMat) <- c("MIBEN", "MIBL", "MLR")
    rownames(outMat) <- c("MSE", "var(MSE)")
    outMat
}# END testFun()

library(parallel)
nReps <- 10
mseList <- mclapply(c(1 : nReps),
                    FUN = testFun,
                    parms = parms,
                    mc.cores = 4)

mseList

?predict.lm


testFun(2, parms)
?predict.lm

dens0 <- density(dat3$y)
dens1 <- density(rowMeans(predOut1))
dens2 <- density(rowMeans(predOut2))
dens3 <- density(predOut3)

yLim <- c(0, max(dens0$y, dens1$y, dens2$y, dens3$y))

plot(dens0, ylim = yLim)
lines(dens1, col = "red")
lines(dens2, col = "blue")
lines(dens3, col = "green")

ls(testOut)

plot(testOut2$lambdaHistory[[1]][ , 1], type = "l")









mod1 <- paste(
    paste0("F",
           colnames(dat1),
           " =~ 1*",
           colnames(dat1),
           "\n"),
    collapse = "")

cat(mod1)

## Estimate the sufficient statistics with FIML:
out1 <- lavaan(model           = mod1,
               data            = dat1,
               int.ov.free     = FALSE,
               int.lv.free     = TRUE,
               auto.var        = TRUE,
               auto.fix.single = TRUE,
               missing         = "fiml")
        
## Store the item means and scales:
dataMeans  <- as.vector(inspect(out1, "coef")$alpha)
dataScales <- sqrt(diag(inspect(out1, "coef")$psi))

(dataMeans - colMeans(dat1)) / colMeans(dat1)
(dataScales - apply(dat1, 2, sd)) / apply(dat1, 2, sd)

names(env$dataMeans) <- names(env$dataScales) <- colnames(env$rawData)
