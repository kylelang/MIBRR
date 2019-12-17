### Title:    Run Tests to Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2019-01-01
### Modified: 2019-12-09

rm(list = ls(all = TRUE))
    
library(lars)
library(MIBRR)
library(ggplot2)

## Prep the data:
data(diabetes)

y    <- diabetes$y
X    <- diabetes$x
dat1 <- data.frame(cbind(y, X))

## Run the BEN models:
out1 <- ben(data        = dat1,
            y           = "y",
            X           = colnames(X),
            iterations  = c(250, 100, 50),
            sampleSizes = list(rep(100, 2), rep(5000, 2), rep(5000, 2)),
            nChains     = 4,
            control     = list(
                lambda1Starts = 1.0 + runif(1, -0.25, 0.25),
                lambda2Starts = 15.0 + runif(1, -5, 5),
                optMethod     = "L-BFGS-B"
            )
            )

## Plot MCEM chains:
MIBRR:::plotLambda(out1, "y")

## Plot loglikelihood trace:
MIBRR:::plotLambda(out1, "y", TRUE)

## Choose rep:
rp <- 1

## Extract parameters:
pars <- getParams(out1, "y")
med  <- apply(pars$beta, 2, median)
ci   <- apply(pars$beta, 2, quantile, probs = c(0.025, 0.975))

## Extract lambda estimate:
pars$lambda1
pars$lambda2

## Get the least squares fit:
#dat2 <- data.frame(cbind(y, scale(X)))
fit  <- lm(y ~ ., data = dat1)
summary(fit)

## Create forest plot of estimates:
dat2 <- data.frame(var = names(med),
                   ben = med,
                   lm  = coef(fit),
                   ciL = ci[1, ],
                   ciU = ci[2, ])
dat2 <- dat2[-1, ]

p1 <- ggplot(data = dat2, mapping = aes(x = ben, y = var)) +
    geom_point() +
    theme_classic() +
    geom_vline(xintercept = 0, linetype = "dotted") +
    scale_y_discrete(limits = rev(dat2$var))
p2 <- p1 + geom_segment(mapping = aes(x = ciL, xend = ciU, y = var, yend = var))
p2 + geom_point(mapping = aes(x = lm, y = var), color = "red")

## Compute (relative) L1 norm of betas:
lmL1 <- sum(abs(coef(fit)[-1]))
benL1 <- sum(abs(med[-1]))

lmL1
benL1

benL1 / lmL1


###--------------------------------------------------------------------------###


## Run fully Bayesian BL models:
out2 <- list()
for(rp in 1 : 4) 
    out2[[rp]] <- bl(data         = dat1,
                     y            = "y",
                     X            = colnames(X),
                     doMcem       = FALSE,
                     sampleSizes  = c(1000, 1000),
                     lam1PriorPar = c(1.0, 0.0041)
                     )
## Choose rep:
rp <- 1

## Extract parameters:
pars <- getParams(out2[[rp]], "y")
med  <- apply(pars$beta, 2, median)
ci   <- apply(pars$beta, 2, quantile, probs = c(0.025, 0.975))

## Summarize the posterior distribution of lambda:
lamMed <- median(pars$lambda[1, ])
lamCi  <- quantile(pars$lambda[1, ], probs = c(0.025, 0.975))

lamMed
lamCi

## Create forest plot of estimates:
dat2 <- data.frame(var = names(med),
                   bl  = med,
                   lm  = coef(fit),
                   ciL = ci[1, ],
                   ciU = ci[2, ])
dat2 <- dat2[-1, ]

p1 <- ggplot(data = dat2, mapping = aes(x = bl, y = var)) +
    geom_point() +
    theme_classic() +
    geom_vline(xintercept = 0, linetype = "dotted") +
    scale_y_discrete(limits = rev(dat2$var))
p2 <- p1 + geom_segment(mapping = aes(x = ciL, xend = ciU, y = var, yend = var))
p2 + geom_point(mapping = aes(x = lm, y = var), color = "red")

## Compute (relative) L1 norm of betas:
lmL1 <- sum(abs(coef(fit)[-1]))
blL1 <- sum(abs(med[-1]))

lmL1
blL1

blL1 / lmL1

?bl

