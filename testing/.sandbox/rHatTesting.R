

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
    
    subChainLen <- floor(length(sims[[1]]) / 2)
    nSubChains  <- length(sims) * 2

    subChainLen
    nSubChains

    do.call(cbind, lapply(sims, matrix, nrow = subChainLen))

    extra <- length(sims[[1]]) %% 2

    
    
    if(length(sims) %% nSubChains == 0) {
        simsMat <- matrix(sims, ncol = nSubChains)
    } else {
        simsMat <- matrix(
            sims[1 : (length(sims) - (nSubChains - 1))],
            ncol = nSubChains
        )
    }

    length(simsMat)
    
    wMean     <- colMeans(simsMat)
    grandMean <- mean(simsMat)
    
    wVar <- mean(apply(simsMat, 2, var))
    bVar <- (subChainLen / (nSubChains - 1)) *
        sum((wMean - grandMean)^2)
    tVar <- ((subChainLen - 1) / subChainLen) * wVar +
        (1 / subChainLen) * bVar
    
    sqrt(tVar / wVar)
},
