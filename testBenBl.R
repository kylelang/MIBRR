### Title:    Build R Packages
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2016-NOV-05
### Purpose:  Script to help build R packages

rm(list = ls(all = TRUE))

install.packages("glmnet", repos = "http://cloud.r-project.org")

library(RcppEigen)

system("rm source/mibrr/src/RcppExports.cpp \
        rm source/mibrr/R/RcppExports.R \
        rm source/mibrr/src/*.o source/mibrr/src/*.so")
Rcpp::compileAttributes("source/mibrr")
install.packages("source/mibrr", repos = NULL, type = "source")

library(mibrr)
library(mitools)
library(glmnet)
library(parallel)

## Test BEN and BL:
nReps          <- 500
alpha          <- 50
nPreds         <- 10

parms                 <- list()
parms$nObs            <- 1000
parms$nPreds          <- nPreds
parms$r2              <- 0.5
parms$collin          <- 0.3
parms$beta            <- matrix(
    c(alpha, runif(ceiling(nPreds / 2), -5, 5), rep(0, floor(nPreds / 2)))
)
parms$means           <- runif(nPreds, -10, 10)
parms$scales          <- runif(nPreds, 1, 50)
parms$simpleInt       <- FALSE
parms$verbose         <- FALSE
parms$adaptScales     <- TRUE
parms$iterations      <- c(500, 25)
parms$latentStructure <- FALSE
parms$itemsPerFactor  <- 4
parms$itemReliability <- 0.8
parms$twoPhaseOpt     <- TRUE
parms$center          <- TRUE
parms$scale           <- TRUE


testFun <- function(rp, parms) {
    nObs            <- parms$nObs
    nPreds          <- parms$nPreds
    r2              <- parms$r2
    collin          <- parms$collin
    beta            <- parms$beta
    means           <- parms$means
    scales          <- parms$scales
    latentStructure <- parms$latentStructure
    itemsPerFactor  <- parms$itemsPerFactor
    itemReliability <- parms$itemReliability

    outList <- list()
    outList$data <- simulateData(nObs            = nObs,
                                 nPreds          = nPreds,
                                 r2              = r2,
                                 collin          = collin,
                                 beta            = beta,
                                 means           = means,
                                 scales          = scales,
                                 latentStructure = latentStructure,
                                 itemsPerFactor  = itemsPerFactor,
                                 itemReliability = itemReliability)
    
    nIVs <- ifelse(latentStructure, nPreds * itemsPerFactor, nPreds)

    outList$ben <- ben(data       = outList$data,
                       y          = "y",
                       X          = paste0("x", c(1 : nIVs)),
                       iterations = parms$iterations,
                       verbose    = parms$verbose,
                       control    =
                           list(center          = parms$center,
                                scale           = parms$scale,
                                adaptScales     = parms$adaptScales,
                                simpleIntercept = parms$simpleInt,
                                twoPhaseOpt     = parms$twoPhaseOpt)
                       )
    
    outList$bl <- bl(data       = outList$data,
                     y          = "y",
                     X          = paste0("x", c(1 : nIVs)),
                     iterations = parms$iterations,
                     verbose    = parms$verbose,
                     control    =
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


outList <- mclapply(X        = c(1 : nReps),
                    FUN      = testFun,
                    parms    = parms,
                    mc.cores = 4)

getwd()

saveRDS(outList, "../../predTestOutList.rds")

















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
