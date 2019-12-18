
rm(list = ls(all = TRUE))

                                        #library(coda)
library(MIBRR)
library(lars)
library(mvtnorm)

data(mibrrExampleData)
dim(mibrrExampleData)

sigma       <- matrix(0.3, 30, 30)
diag(sigma) <- 1.0

dat1           <- as.data.frame(rmvnorm(20, rep(0, 30), sigma))
colnames(dat1) <- paste0("x", 1 : 30)

## MCEM estimation:
blOut <- bl(data    = dat1,
            y       = "x1",
            X       = paste0("x", 2 : 30),
            nChains = 2)

benOut <- ben(data    = dat1,
              y       = "x1",
              X       = paste0("x", 2 : 30),
              nChains = 2)

plotLambda(benOut, "x1")

lams <- getLambda(benOut, "x1")

l1 <- lapply(lams, function(x, n) tail(x$lambda1, n), n = 100)
l2 <- lapply(lams, function(x, n) tail(x$lambda2, n), n = 100)
ll <- lapply(lams, function(x, n) tail(x$logLik, n), n = 100)

l1 <- unlist(lapply(l1, MIBRR:::split), recursive = FALSE)
l2 <- unlist(lapply(l2, MIBRR:::split), recursive = FALSE)
ll <- unlist(lapply(ll, MIBRR:::split), recursive = FALSE)

mct <- function(sams) {
    ## Compute the medians and ranges of each sample:
    med <- sapply(sams, median)
    ran <- sapply(sams, range)
    
    ## Test that each median is within the range of every sample:
    all(
        sapply(X   = med,
               FUN = function(m, r) m >= r[1, ] & m <= r[2, ],
               r   = ran)
    )
}

unlist(lapply(benOut$getSamples("logLik", "x1"), MIBRR:::split), recursive = FALSE)

mct(l1)
mct(l2)
mct(ll)

all(checks)

unlist(l1, FALSE)

l1[[1]]

l1[[1]]

l1
l1
ll

mibenOut

mibenOut$rHats

## Prep the data:
data(diabetes)

y    <- diabetes$y
X    <- diabetes$x
dat1 <- data.frame(cbind(y, X))

## Run the BEN models:
out1 <- ben(data        = dat1,
            y           = "y",
            X           = colnames(X),
            iterations  = c(500, 500),
            sampleSizes = list(rep(250, 2), rep(500, 2), rep(1000, 2)),
            nChains     = 4,
            control     = list(
                lambda1Starts = 1.0 + runif(1, -0.25, 0.25),
                lambda2Starts = 15.0 + runif(1, -5, 5),
                optMethod     = "L-BFGS-B"
            )
            )

rHats <- getField(out1, "rHats")
rHats

pars <- getParams(out1, "y", FALSE)

names(pars)

beta <- pars$beta

tmp <- beta[[1]]
x <- pars$sigma[[1]][1 : 499]

split <- function(x) {
    ## Make sure we're working with a matrix:
    if(!is.matrix(x)) x <- as.matrix(x)

    ## Make sure N is a multiple of two:
    extra <- nrow(x) %% 2
    if(extra != 0) x <- as.matrix(x[(extra + 1) : nrow(x), ])

    ## Define subchain length:
    n <- nrow(x) / 2

    ## Split the sample:
    list(x[1 : n, ],
         x[(n + 1) : nrow(x), ]
         )
}

prepSams <- function(sams, split = TRUE) {
    ## Split each sample into two subchain samples:
    if(split)
        sams <- do.call(c, lapply(sams, split))
    
    ## Convert samples into an mcmc.list object:
    mcmc.list(lapply(sams, mcmc))
}

tmp <- split(pars$sigma[[1]])
tmp

tmp <- prepSam(pars$lambda1, FALSE)
tmp

n <- out1$iterations[2]

iterations <- out1$iterations
iterations[2] <- 50

prepMcemChains <- function(sams) {
    sams <- lapply(X   = sams,
                   FUN = function(x, n) x[(length(x) - n + 1) : length(x)],
                   n   = ceiling(iterations[2] / 2)
                   )
    
    prepSams(sams)
}

tmp <- prepMcemChains(pars$lambda1)

tmp

m <- mean(unlist(tmp))
s <- sd(unlist(tmp))

lapply(tmp, function(x, m, s) (x - m) / s, m = m, s = s)

rHats <- gelman.diag(tmp, transform = TRUE)
rHats

geweke.diag(tmp)

## Plot MCEM chains:
par(mfrow = c(1, 2))

lam1 <- getParams(out1, "y", mix = FALSE)$lambda1
lam2 <- getParams(out1, "y", mix = FALSE)$lambda2

cols <- rainbow(length(lam1))

plot(lam1[[1]], type = "l", col = cols[1])
for(i in 2 : length(lam1))
    lines(lam1[[i]], col = cols[i])

plot(lam2[[1]], type = "l", col = cols[1])
for(i in 2 : length(lam2))
    lines(lam2[[i]], col = cols[i])

length(tmp)

par(mfrow = c(2, 4))
lapply(tmp, function(x) plot(density(log(x))))

class(rHats)
rHats$psrf[ , 1]
rHats$psrf[ , 2]

geweke.diag(tmp)
heidel.diag(tmp)

length(tmp)

dim(tmp[[1]])

head(tmp[[1]])

length(beta)

tmp <- mcmc(beta[[1]])

tmp

mcmcBeta <- mcmc.list(lapply(beta, mcmc))

tmp <- gelman.diag(mcmcBeta)

ls(tmp)


sims <- rep(1 : 3, each = 19)
sims <- rep(1, 19)

nChains <- 3

sims

n <- 9

sims <- list(1 : n,
             1 : n,
             1 : n)

sims

calcRHat <- function(sims) {
    "Compute a single split-chain Potential Scale Reduction Factor"
    if(!is.list(sims)) stop("'sims' must be a list")

    ## Compute the number of sub-chains:
    nSubChains  <- length(sims) * 2

    ## Make sure the number of samples is a multiple of two:
    extra <- length(sims[[1]]) %% 2
    sims  <- lapply(sims, function(x, y) x[(y + 1) : length(x)], y = extra)

    ## Define the length of each sub-chain:
    subChainLen <- length(sims[[1]]) / 2

    ## Restructure the chains to facilitate variance calculations:
    simsMat <- do.call(cbind, lapply(sims, matrix, nrow = subChainLen))
    
    wMean     <- colMeans(simsMat)
    grandMean <- mean(simsMat)
    
    wVar <- mean(apply(simsMat, 2, var))
    bVar <- (subChainLen / (nSubChains - 1)) * sum((wMean - grandMean)^2)
    tVar <- ((subChainLen - 1) / subChainLen) * wVar + (1 / subChainLen) * bVar
    
    sqrt(tVar / wVar)

    w <- mean(apply(simsMat, 2, var))
    b <- mean(apply(simsMat, 1, var))

    w
    b

    b / w
},
