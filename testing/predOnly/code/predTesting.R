### Title:    Test MIBEN/MIBL Prediction
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2016-MAY-15

rm(list = ls(all = TRUE))

outDir <- "../output/"

library(mice)
library(mitools)
library(mvtnorm)
library(mibrr)
library(parallel)
library(Matrix)

## Test BEN and BL:
alpha  <- 50
nPreds <- 90
nReps  <- 500

condMat <- expand.grid(center = c(TRUE, FALSE), scale = c(TRUE, FALSE))

means  <- runif(nPreds, 0, 100)
scales <- runif(nPreds, 1, 100)
collin <- 0.3

outMeans <- matrix(NA, 4, 3)
outSds <- matrix(NA, 4, 3)
outMse <- list()
for(i in 1 : nrow(condMat)) {
    print(paste0("Doing condition set ", i))
    print(" ")
    
    parms <- list()
    parms$nObs        <- 100
    parms$nPreds      <- nPreds
    parms$r2          <- 0.5
    parms$collin      <- collin
                                        #parms$beta        <- matrix(
                                        #    c(alpha, runif(nPreds / 2, 0.3, 0.6), rep(0, nPreds / 2))
                                        #)
    parms$beta        <- matrix(c(alpha, runif(nPreds, 0.3, 0.6)))
    parms$means       <- means
    parms$scales      <- scales
    parms$center      <- condMat[i, "center"]
    parms$scale       <- condMat[i, "scale"]
    parms$useClassic  <- FALSE
    parms$simpleInt   <- FALSE
    parms$verbose     <- FALSE
    parms$postN       <- 1000
    parms$approxIters <- 200
    parms$approxN     <- 100
    parms$tuneIters   <- 20
    parms$tuneN       <- 250

    mseList <- mclapply(c(1 : nReps),
                        FUN = testFun,
                        parms = parms,
                        mc.cores = 4)
  
    mseMat <- do.call(rbind, mseList)
    outMse[[i]] <- mseMat
    outMeans[i, ] <- colMeans(mseMat)
    outSds[i, ] <- apply(mseMat, 2, sd)
}


outList <- list(
    mean = cbind(condMat, outMeans),
    sd =  cbind(condMat, outSds)
)

outList

outList$mean[ , c(3, 4)] / outList$mean[ , 5]
outList$sd[ , c(3, 4)] / outList$sd[ , 5]

saveRDS(outList,
        file = paste0(outDir, "dense_90pred_outList1.rds")
        )

##> cbind(condMat, outMeans)
##  center scale        1        2        3
##1   TRUE  TRUE 30.30673 31.31431 23.86253
##2  FALSE  TRUE 28.82615 29.18330 22.49048
##3   TRUE FALSE 39.95802 39.04987 30.71749
##4  FALSE FALSE 29.02999 28.63054 22.85032
##> cbind(condMat, outSds)
##  center scale        1        2        3
##1   TRUE  TRUE 12.10952 15.52743 4.647798
##2  FALSE  TRUE 11.49169 12.17685 4.394110
##3   TRUE FALSE 17.08081 18.00464 6.517262
##4  FALSE FALSE 11.52794 11.17665 4.767104

testFun <- function(rp, parms) {
    print(paste0("Doing replication ", rp))
    
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
    
    testOut <- ben(rawData         = dat1,
                   y               = "y",
                   X               = paste0("x", c(1 : nPreds)),
                   mcemApproxIters = parms$approxIters,
                   mcemApproxN     = parms$approxN,
                   mcemTuneIters   = parms$tuneIters,
                   mcemTuneN       = parms$tuneN,
                   mcemPostN       = parms$postN,
                   verbose         = parms$verbose,
                   control         =
                       list(center          = parms$center,
                            scale           = parms$scale,
                            adaptScales     = FALSE,
                            regIntercept    = FALSE,
                            useClassic      = parms$useClassic,
                            simpleIntercept = parms$simpleInt)
                   )
    
    testOut2 <- bl(rawData         = dat1,
                   y               = "y",
                   X               = paste0("x", c(1 : nPreds)),
                   mcemApproxIters = parms$approxIters,
                   mcemApproxN     = parms$approxN,
                   mcemTuneIters   = parms$tuneIters,
                   mcemTuneN       = parms$tuneN,
                   mcemPostN       = parms$postN,
                   verbose         = parms$verbose,
                   control         =
                       list(center          = parms$center,
                            scale           = parms$scale,
                            adaptScales     = FALSE,
                            regIntercept    = FALSE,
                            useClassic      = parms$useClassic,
                            simpleIntercept = parms$simpleInt)
                   )

    testForm <-
        as.formula(paste0("y ~ ", paste0("x", c(1 : nPreds), collapse = " + ")))
    lmOut <- lm(testForm, data = dat1)
    
    testDat <- simulateData(nObs   = nObs,
                            nPreds = nPreds,
                            r2     = r2,
                            collin = collin,
                            beta   = beta,
                            means  = means,
                            scales = scales)

    if(parms$center) {
        testDat2 <- data.frame(scale(testDat, scale = FALSE))
        extra <- mean(dat1$y)
    } else {
        testDat2 <- testDat
        extra <- 0.0
    }
    if(parms$scale) {
        tmpSds    <- apply(testDat2, 2, sd)
        tmpSds[1] <- 1.0
        testDat2  <- data.frame(
            scale(testDat2, center = FALSE, scale = tmpSds)
        )
    }
    
    predOut1 <- predictMibrr(object = testOut,
                             newData = as.matrix(testDat2[ , -1])
                             ) + extra
    predOut2 <- predictMibrr(object = testOut2,
                             newData = as.matrix(testDat2[ , -1])
                             ) + extra
    predOut3 <- predict(lmOut, newdata = testDat)
    
    mseVec <- c(
        mean((predOut1 - dat1$y)^2),
        mean((predOut2 - dat1$y)^2),
        mean((predOut3 - dat1$y)^2)
    )
    mseVec
}# END testFun()

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
dat1.2 <- data.frame(scale(dat1))
dat1.3 <- data.frame(scale(dat1, scale = FALSE))
dat1.4 <- data.frame(scale(dat1, center = FALSE, scale = apply(dat1, 2, sd)))
dat1.5 <- data.frame(dat1.3[ , 1], dat1.2[ , -1])
dat1.6 <- data.frame(dat1[ , 1], dat1.4[ , -1])

colnames(dat1.2) <- colnames(dat1.3) <- colnames(dat1.4) <- colnames(dat1.5) <-
    colnames(dat1.6) <- colnames(dat1)

dat2 <- simulateData(nObs   = nObs,
                     nPreds = nPreds,
                     r2     = r2,
                     collin = collin,
                     beta   = beta,
                     means  = means,
                     scales = scales)
dat2.2 <- data.frame(scale(dat2))
dat2.3 <- data.frame(scale(dat2, scale = FALSE))
dat2.4 <- data.frame(scale(dat2, center = FALSE, scale = apply(dat2, 2, sd)))
dat2.5 <- data.frame(dat2.3[ , 1], dat2.2[ , -1])
dat2.6 <- data.frame(dat2[ , 1], dat2.4[ , -1])

colnames(dat2.2) <- colnames(dat2.3) <- colnames(dat2.4) <- colnames(dat2.5) <-
    colnames(dat2.6) <- colnames(dat2)

form1 <- as.formula(paste0("y ~ ", paste0("x", c(1 : 10), collapse = " + ")))

out1 <- lm(form1, data = dat1)   # Raw
out2 <- lm(form1, data = dat1.2) # Standardized
out3 <- lm(form1, data = dat1.3) # Mean Centered
out4 <- lm(form1, data = dat1.4) # Scaled
out5 <- lm(form1, data = dat1.5) # Mix
out6 <- lm(form1, data = dat1.6)

pred1 <- predict(out1, newdata = dat2)
pred2 <- predict(out6, newdata = dat2)
pred3 <- predict(out6, newdata = dat2.6)

dens1 <- density(pred1)
dens2 <- density(pred2)
dens3 <- density(pred3)

yLim <- range(dens1$y, dens2$y, dens3$y)
xLim <- range(dens1$x, dens2$x, dens3$x)

plot(density(pred1), col = "red", ylim = yLim, xlim = xLim)
lines(density(pred2), col = "blue")
lines(density(pred3), col = "green")
