### Title:    Test MIBRR Package
### Author:   Kyle M. Lang
### Created:  2014-DEC-07
### Modified: 2017-NOV-27

rm(list = ls(all = TRUE))

                                        #install.packages("HyperbolicDist", repos = "http://cloud.r-project.org")

library(mitools)
library(psych)
library(MIBRR)
library(devtools)
library(parallel)
library(MCMCpack)
library(statmod)
library(HyperbolicDist)

install_github("kylelang/SURF/source/SURF")
library(SURF)

saveDate <- format(Sys.time(), "%Y%m%d")

###---------------------------------------------------------------------------###

### Prepare Data for Testing ###

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

targets  <- list(mar = paste0(c("A", "E"), rep(c(1 : 5), each = 2)),
                 mcar = NA,
                 mnar = NA)
pm       <- list(mar = 0.3)
snr      <- list(mar = 5)
marPreds <- c("age",
              "male",
              "finish_hs",
              "some_college",
              "college_grad",
              "graduate_degree")
cn       <- c(targets$mar, marPreds)

dat1 <- bfi2[sample(c(1 : nrow(bfi2)), 500), cn]
dat2 <- imposeMissData(data    = dat1,
                       targets = targets,
                       preds   = marPreds,
                       pm      = pm,
                       snr     = snr)$data 

###---------------------------------------------------------------------------###

### Check fully Bayesian estimation of MIBEN and MIBL ###

mibenOut <- miben(data         = dat2,
                  targetVars   = targets$mar,
                  ignoreVars   = NULL,
                  iterations   = c(50, 10),
                  sampleSizes  = c(1000, 1000),
                  doMcem       = FALSE,
                  lam1PriorPar = c(1.0, 10.0),
                  lam2PriorPar = c(2.0, 0.01),
                  verbose      = TRUE)

mibenImps <- MIBRR::complete(mibenOut, 100)

mibenImps[[1]]

miblOut <- mibl(data       = dat2,
                targetVars = targets$mar,
                ignoreVars = NULL,
                iterations = c(50, 10),
                verbose    = FALSE)

miblImps <- MIBRR::complete(miblOut, 100)


###---------------------------------------------------------------------------###

### Test Random Variate Samplers ###

## MVN sampler:
nObs <- 500000
mvnMu <- rep(10, 3)
mvnSigma <- matrix(5, 3, 3)
diag(mvnSigma) <- 20

out1.1 <- MIBRR:::drawMvn(nObs, mvnMu, mvnSigma)
out1.2 <- rmvnorm(nObs, mvnMu, mvnSigma)

par(mfrow = c(1, 3))
for(i in 1 : length(mvnMu)) {
    plot(density(out1.1[ , i]), col = "red")
    lines(density(out1.2[ , i]), col = "blue")
}

## Inverse gamma sampler:
gamShape <- 10
gamScale <- 10

out2.1 <- MIBRR:::drawInvGamma(nObs, gamShape, gamScale)
out2.2 <- rinvgamma(nObs, gamShape, gamScale)

plot(density(out2.1), col = "red")
lines(density(out2.2), col = "blue")

## Inverse gaussian sampler:
igMu  <- 1
igLam <- 2

out3.1 <- MIBRR:::drawInvGauss(nObs, igMu, igLam)
out3.2 <- rinvgauss(nObs, igMu, igLam)

plot(density(out3.1), col = "red")
lines(density(out3.2), col = "blue")

## GIG sampler:
gigLam <- 1
gigChi <- 2
gigPsi <- 2

out4.1 <- MIBRR:::drawGig(nObs, gigLam, gigChi, gigPsi)
out4.2 <- rgig(nObs, c(gigLam, gigChi, gigPsi))

plot(density(out4.1), col = "red")
lines(density(out4.2), col = "blue")


## Incomplte gamma calculation:
incGamShape <- 10
incGamCut   <- 5

out4.1 <- MIBRR:::calcIncGamma(incGamShape, incGamCut, FALSE)
out4.2 <- pgamma(q     = incGamCut,
                 shape = incGamShape,
                 lower = FALSE) * gamma(incGamShape)

out4.1 - out4.2
