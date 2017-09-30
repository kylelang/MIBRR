### Title:    Test BEN and BL
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2016-NOV-11
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))

                                        #install.packages("glmnet", repos = "http://cloud.r-project.org")

                                        #library(RcppEigen)

                                        #softDir <- "~/data/software/development/mibrrProject/active/mibrr/"

                                        #system(
                                        #    paste0("rm ", softDir, "source/mibrr/src/RcppExports.cpp\n",
                                        #           "rm ", softDir, "source/mibrr/R/RcppExports.R\n",
                                        #           "rm ", softDir, "source/mibrr/src/*.o\n",
                                        #           "rm ", softDir, "source/mibrr/src/*.so\n")
                                        #)
                                        #Rcpp::compileAttributes(paste0(softDir, "source/mibrr"))
                                        #install.packages(paste0(softDir, "source/mibrr"),
                                        #                 repos = NULL,
                                        #                 type  = "source")

library(mibrr)
library(mitools)
library(glmnet)
library(parallel)

## Test BEN and BL:
outDir <- "../output/"
alpha  <- 50
nPreds <- 40

parms                 <- list()
parms$nObs            <- 50
parms$nPreds          <- nPreds
parms$r2              <- 0.2
parms$collin          <- 0.5
parms$beta            <- matrix(
                                        #c(alpha, runif(ceiling(nPreds / 2), 5, 10), rep(0, floor(nPreds / 2)))
                                        #c(alpha, runif(nPreds, 5, 10))
    c(alpha, rep(0.3, nPreds))
)
parms$means           <- rep(0, nPreds) #runif(nPreds, -10, 10)
parms$scales          <- rep(1, nPreds) #runif(nPreds, 1, 50)
parms$simpleInt       <- FALSE
parms$verbose         <- FALSE
parms$adaptScales     <- TRUE
parms$iterations      <- c(1000, 100)
parms$sampleSizes     <- c(50, 100, 1000)
parms$latentStructure <- FALSE
parms$itemsPerFactor  <- 4
                                        #parms$itemReliability <- 1.0
parms$twoPhaseOpt     <- TRUE
parms$center          <- TRUE
parms$scale           <- TRUE
parms$relRange        <- c(0.5, 0.75)

simulateData2 <- function(nObs,
                          nPreds,
                          r2,
                          collin,
                          beta,
                          means           = 0,
                          scales          = 1,
                          latentStructure = FALSE,
                          itemsPerFactor  = 1,
                                        #itemReliability = NULL,
                          relRange        = c(0.5, 0.75))
{
    if(length(means) == 1) means <- rep(means, nPreds)
    
    w1 <- matrix(scales, nPreds, nPreds)
    w2 <- matrix(scales, nPreds, nPreds, byrow = TRUE)
    
    maxCov <- w1*w2
    
    sigma <- maxCov * collin
    diag(sigma) <- scales^2
    
    X <- cbind(1, rmvnorm(nObs, means, sigma))

    eta <- X %*% beta
    sigmaY <- (var(eta) / r2) - var(eta)
    y <- eta + rnorm(nObs, 0, sqrt(sigmaY))

    if(latentStructure) {
        nItems   <- nPreds * itemsPerFactor
        loadings <- matrix(0, nItems, nPreds)
        
        for(m in 1 : nPreds) {
            for(n in 1 : itemsPerFactor) {
                offset <- (m - 1) * itemsPerFactor
                loadings[n + offset, m] <-
                    sqrt(runif(1, relRange[1], relRange[2]))
            }
        }
                
        theta <- diag(1 - loadings[loadings != 0]^2)
    
        X <- X[ , -1] %*% t(loadings) + rmvnorm(nObs, rep(0, nItems), theta)
        
        outDat <- data.frame(y, X)
        colnames(outDat) <- c("y", paste0("x", c(1 : nItems)))
    } else {
        outDat <- data.frame(y, X[ , -1])
        colnames(outDat) <- c("y", paste0("x", c(1 : nPreds)))
    }
    outDat
}



testFun <- function(rp, parms) {
    print(paste0("Doing replication ", rp))
    
    nObs            <- parms$nObs
    nPreds          <- parms$nPreds
    r2              <- parms$r2
    collin          <- parms$collin
    beta            <- parms$beta
    means           <- parms$means
    scales          <- parms$scales
    latentStructure <- parms$latentStructure
    itemsPerFactor  <- parms$itemsPerFactor
                                        #itemReliability <- parms$itemReliability
    relRange        <- parms$relRange

    outList <- list()
    outList$data <- simulateData2(nObs            = nObs,
                                  nPreds          = nPreds,
                                  r2              = r2,
                                  collin          = collin,
                                  beta            = beta,
                                  means           = means,
                                  scales          = scales,
                                  latentStructure = latentStructure,
                                  itemsPerFactor  = itemsPerFactor,
                                        #itemReliability = itemReliability
                                  relRange        = relRange)
    
    nIVs <- ifelse(latentStructure, nPreds * itemsPerFactor, nPreds)

    outList$ben <- ben(data        = outList$data,
                       y           = "y",
                       X           = paste0("x", c(1 : nIVs)),
                       iterations  = parms$iterations,
                       sampleSizes = parms$sampleSizes,
                       verbose     = parms$verbose,
                       control     =
                           list(center          = parms$center,
                                scale           = parms$scale,
                                adaptScales     = parms$adaptScales,
                                simpleIntercept = parms$simpleInt,
                                twoPhaseOpt     = parms$twoPhaseOpt)
                       )
    
    outList$bl <- bl(data        = outList$data,
                     y           = "y",
                     X           = paste0("x", c(1 : nIVs)),
                     iterations  = parms$iterations,
                     sampleSizes = parms$sampleSizes,
                     verbose     = parms$verbose,
                     control     =
                         list(center          = parms$center,
                              scale           = parms$scale,
                              adaptScales     = parms$adaptScales,
                              simpleIntercept = parms$simpleInt)
                     )
    
    testForm <-
        as.formula(paste0("y ~ ", paste0("x", c(1 : nIVs), collapse = " + ")))

                                        #tmpDat <- as.data.frame(scale(outList$data, scale = FALSE))
    outList$lm <- lm(testForm, data = outList$data)
    outList
}


mseFun <- function(x, nReps) {
    nIVs <- ifelse(parms$latentStructure, parms$nPreds * parms$itemsPerFactor,
                   parms$nPreds)

    outMat <- matrix(NA, nReps, length(outList[[1]]) - 1)
    for(rp in 1 : nReps) {
        testData <- simulateData2(nObs            = parms$nObs,
                                  nPreds          = parms$nPreds,
                                  r2              = parms$r2,
                                  collin          = parms$collin,
                                  beta            = parms$beta,
                                  means           = parms$means,
                                  scales          = parms$scales,
                                  latentStructure = parms$latentStructure,
                                  itemsPerFactor  = parms$itemsPerFactor,
                                        #itemReliability = parms$itemReliability
                                  relRange        = parms$relRange)
        
        newX <- scale(testData[ , paste0("x", c(1 : nIVs))])
        yMean <- mean(testData$y)
        
        nDraws <- length(x$ben$params$y$sigma)
        
        benPPD <- predictMibrr(object    = x$ben,
                               newData   = newX,
                               targetVar = "y",
                               nDraws    = nDraws)
        benPred <- rowMeans(benPPD) + yMean
        
        blPPD <- predictMibrr(object    = x$bl,
                              newData   = newX,
                              targetVar = "y",
                              nDraws    = nDraws)
        blPred <- rowMeans(blPPD) + yMean
        
        lmPred    <- predict(x$lm, newdata = testData[ , -1])
        
        outMat[rp, ] <- c(mean((benPred   - testData$y)^2),
                          mean((blPred    - testData$y)^2),
                          mean((lmPred    - testData$y)^2)
                          )
        colnames(outMat) <- c("MIBEN", "MIBL", "MLR")
    }
    colMeans(outMat)
}


outList <- testFun(1, parms)

tmp1 <- outList$ben$lambdaHistory$y
tmp2 <- outList$bl$lambdaHistory$y

par(mfrow = c(1, 3))
plot(tmp1[ , "lambda1"], type = "l")
plot(tmp1[ , "lambda2"], type = "l")
plot(tmp2, type = "l")

outList <- mclapply(X        = c(1 : 100),
                    FUN      = testFun,
                    parms    = parms,
                    mc.cores = 4)
saveRDS(outList, paste0(outDir, "predTestOutList19.rds"))


                                        #outList <- readRDS(paste0(outDir, "predTestOutList.rds"))


mean(benPred)
mean(blPred)
mean(lmPred)
mean(testData$y)

d1 <- density(testData$y)
d2 <- density(benPred)
d3 <- density(blPred)
d4 <- density(lmPred)

var(benPred)
var(blPred)
var(lmPred)
var(testData$y)

yLim <- range(d1$y, d2$y, d3$y, d4$y)
xLim <- range(d1$x, d2$x, d3$x, d4$x)

par(mfrow = c(1, 1))
plot(d1, ylim = yLim, xlim = xLim)
lines(d2, col = "red")
lines(d3, col = "blue")
lines(d4, col = "green")

mseFrame <- do.call(rbind, lapply(outList, mseFun, nReps = 50))

colMeans(mseFrame)
apply(mseFrame, 2, median)
apply(mseFrame, 2, var)

cor(benPred, testData$y)
cor(blPred, testData$y)
cor(lmPred, testData$y)

plot(testData$y, benPred)
plot(testData$y, blPred)
plot(testData$y, lmPred)




out1 <- testFun(1, parms)

benOut <- out1$ben
blOut  <- out1$bl
lmOut  <- out1$lm

pars1 <- benOut$params$y

colMeans(pars1$beta)
colMeans(pars1$tau)
mean(pars1$sigma)


pars2 <- blOut$params$y

colMeans(pars2$beta)
colMeans(pars2$tau)
mean(pars2$sigma)

dat1 <- out1$data

benPPD <- predictMibrr(object    = benOut,
                       newData   = scale(dat1[ , -1]),
                       targetVar = "y",
                       nDraws    = nrow(benOut$params$y$beta)
                       ) + mean(dat1$y)

blPPD <- predictMibrr(object    = blOut,
                      newData   = scale(dat1[ , -1]),
                      targetVar = "y",
                      nDraws    = nrow(benOut$params$y$beta)
                      ) + mean(dat1$y)

lmPPD <- predict(lmOut, newData = dat1[ , -1])

benDens <- density(rowMeans(benPPD))
blDens  <- density(rowMeans(blPPD))
lmDens  <- density(lmPPD)
yDens   <- density(scale(dat1$y, scale = FALSE))

yLim <- range(benDens$y, blDens$y, lmDens$y, yDens$y)
xLim <- range(benDens$x, blDens$x, lmDens$x, yDens$x)

par(mfrow = c(1, 1))
plot(benDens, col = "red", ylim = yLim, xlim = xLim)
lines(blDens, col = "blue")
lines(lmDens, col = "green")
lines(yDens)

mean(benPPD)
mean(dat1$y)
mean(lmPPD)
mean(blPPD)

plot(enOut)

?plot.glmnet

?cv.glmnet
?glmnet
?enet

x2 <- apply(dat1, 2, as.numeric)

plot(lassoOut)
?plot.glmnet

lassoOut <- cv.glmnet(x = as.matrix(dat1[ , paste0("x", c(1 : nIVs))]),
                      y = dat1$y,
                      foldid = c(1 : nrow(dat1)),
                      alpha = 1)

ridgeOut <- cv.glmnet(x = as.matrix(dat1[ , paste0("x", c(1 : nIVs))]),
                      y = dat1$y,
                      foldid = c(1 : nrow(dat1)),
                      alpha = 0)

X <- as.matrix(dat1[ , paste0("x", c(1 : nIVs))])
y <- dat1$y


(is.na(cvmVec))
lamVec

warnings()

ls(ridgeOut)

ridgeOut$glmnet.fit
ridgeOut$lambda.min
lassoOut$lambda.min

library(enet)

out1 <- testFun(1, parms)
out1

mean(varPred1)
mean(varPred2)

library(parallel)
nReps <- 10
mseList <- mclapply(c(1 : nReps),
                    FUN = testFun,
                    parms = parms,
                    mc.cores = 4)

tmp1 <- benOut$lambdaHistory$y
tmp2 <- blOut$lambdaHistory$y

par(mfrow = c(1, 3))
plot(tmp1[ , "lambda1"], type = "l")
plot(tmp1[ , "lambda2"], type = "l")
plot(tmp2, type = "l")


nObs <- 100000

x <- rnorm(nObs, 0, 10)

alpha <- 100
beta <- 50
r2 <- 0.5

eta <- alpha + beta*x
sigma <- (var(eta) / r2) - var(eta)

y <- eta + rnorm(nObs, 0, sqrt(sigma))

summary(lm(y ~ x))

mean(y)
mean(x)

y2 <- y - 200

out1 <- lm(y ~ x)
out2 <- lm(y2 ~ x)
out3 <- lm(y ~ x - 1)
out4 <- lm(y2 ~ x - 1)

pred1 <- predict(out1, newdata = data.frame(x))
pred2 <- predict(out2, newdata = data.frame(x))
pred3 <- predict(out3, newdata = data.frame(x))
pred4 <- predict(out4, newdata = data.frame(x))

d0 <- density(y)
d00 <- density(y2)
d1 <- density(pred1)
d2 <- density(pred2)
d3 <- density(pred3)
d4 <- density(pred4)

yLim <- range(d1$y, d2$y, d3$y, d4$y)
xLim <- range(d1$x, d2$x, d3$x, d4$x)

plot(d1, ylim = yLim, xlim = xLim)
lines(d2, col = "red")
lines(d3, col = "blue")
lines(d4, col = "green")
lines(d0, col = "purple")
lines(d00, col = "orange")

mean(y)
mean(y2)
mean(pred1)
mean(pred2)
mean(pred3)
mean(pred4)

coef(out1)[1] + coef(out1)[2] * mean(x)
mean(pred1)

mean(pred2 + 200)
mean(pred4 + 200)

y3 <- y
y3[x > quantile(x, 0.5)] <- NA

mean(y3, na.rm = TRUE)

out5 <- lm(y3 ~ x)
out6 <- lm(y3 ~ x - 1)
summary(out1)
summary(out5)
summary(out6)

pred5 <- predict(out5, newdata = data.frame(x))
pred6 <- predict(out6, newdata = data.frame(x))

d5 <- density(pred5)
d6 <- density(pred6)

yLim <- range(d5$y, d6$y)
xLim <- range(d5$x, d6$x)

plot(d5, ylim = yLim, xlim = xLim)
lines(d6, col = "red")

mean(pred5)
mean(pred6)

y4 <- y3 - mean(y)

x2 <- data.frame(x[!is.na(y3)])
x3 <- data.frame(x[is.na(y3)])
colnames(x2) <- colnames(x3) <- "x2"

y4 <- matrix(na.omit(y3))

mean(y4)

out1 <- lm(y ~ x)
out2 <- lm(y3 ~ x)
out3 <- lm(y4 ~ x2)

summary(out1)
summary(out2)
summary(out3)

pred1 <- predict(out1, data.frame(x))
pred2 <- predict(out3, x2)
pred3 <- predict(out3, x3)


pred2.2 <- cbind(1, x2) %*% matrix(coef(out3))
pred3.2 <- cbind(1, x3) %*% matrix(coef(out3))


d1 <- density(pred1)
d2 <- density(pred2)
d3 <- density(pred3)
d4 <- density(c(pred2, pred3))

yLim <- range(d1$y, d2$y, d3$y)
xLim <- range(d1$x, d2$x, d3$x)

plot(d1, ylim = yLim, xlim = xLim)
lines(d2, col = "red")
lines(d3, col = "blue")
lines(d4, col = "green")

?predict.lm

mean(x2)
mean(x3)

mean(pred2)
mean(pred3)


mean(x2)
mean(x3)

out7 <- lm(y4[!is.na(y4)] ~ x[!is.na(y4)])
out8 <- lm(y4[!is.na(y4)] ~ x[!is.na(y4)] - 1)

pred7 <- predict(out7, newdata = data.frame(x3))
pred8 <- predict(out8, newdata = data.frame(x3))

d7 <- density(pred7)
d8 <- density(pred8)

yLim <- range(d7$y, d8$y)
xLim <- range(d7$x, d8$x)

plot(d7, ylim = yLim, xlim = xLim)
lines(d8, col = "red")

mean(pred7)
mean(pred8)

length(x2)
length(x3)

?lm

pred9.2 <- cbind(1, x2) %*% matrix(coef(out7))
pred10.2 <- x2 %*% matrix(coef(out8))

pred9 <- predict(out7, newdata = data.frame(x2))
pred10 <- predict(out8, newdata = data.frame(x2))

all.equal(matrix(pred9), matrix(pred9.2))
all.equal(matrix(pred10), matrix(pred10.2))

summary(out7)
summary(out8)

d7 <- density(pred7)
d8 <- density(pred8)

mean(pred9)
mean(pred10)

mean(pred7)
mean(pred8)

mean(x2)
mean(x3)

summary(out7)
summary(out9)

mean(pred9)


x <- c(1 : 10)

mean(x)

x2 <- 2 * x

mean(x2)
var(x)
var(x)*4
var(x2)

x <- rnorm(100)
x2 <- 4 * x

mean(x)
mean(x2) / mean(x)

var(x2) / var(x)
