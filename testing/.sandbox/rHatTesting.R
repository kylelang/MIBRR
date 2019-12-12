
library(coda)
library(MIBRR)

data(mibrrExampleData)

## MCEM estimation:
mibenOut <- miben(data       = mibrrExampleData,
                  iterations = c(30, 10),
                  targetVars = c("y", paste0("x", c(1 : 3))),
                  ignoreVars = "idNum",
                  nChains    = 3)

pars <- getParams(mibenOut, "y", FALSE)

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

prepSam <- function(sam) {
    ## Split each sample into two subchain samples:
    tmp <- do.call(c, lapply(sam, split))

    ## Convert samples into an mcmc.list object:
    mcmc.list(lapply(tmp, mcmc))
}

tmp <- split(pars$sigma[[1]])
tmp

tmp <- prepSam(pars$beta)

gelman.diag(tmp)
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
