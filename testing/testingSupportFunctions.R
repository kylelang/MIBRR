### Title:    Testing Support Functions for the MIBRR Package
### Author:   Kyle M. Lang
### Created:  2017-OCT-25
### Modified: 2017-OCT-25


makeRVec <- function(linPred, pm, snr, pattern) {
    ## Generate a linear predictor of missingness:
    noise   <- sd(linPred) * (1 / snr)
    linPred <- linPred + rnorm(length(linPred), 0, noise)

    ## Use a probit model to simulate response propensities:
    probs <- pnorm(scale(linPred))

    ## 'pattern' selects which part of the target's distribution is missing:
    if(pattern == "random")
        pattern <- sample(c("low", "high", "tails", "center"), size = 1)

    ## Generate the missignness indicator:
    rVec <- switch(pattern,
                   low    = probs <= pm,
                   high   = probs >= (1 - pm),
                   tails  = probs <= (pm / 2) | probs >= (1 - (pm / 2)),
                   center = probs >= (0.5 - (pm / 2)) &
                       probs <= (0.5 + (pm / 2)),
                   stop("Please provide a valid 'pattern'")
                   )
    list(rVec = rVec, pattern = pattern)
}# END makeRVec()

        

imposeMissing <- function(data, targets, preds, pm, snr, pattern = "random") {
    ## Which mechanisms should be simulated?
    mechFlag <- !is.na(targets)

    if(length(pattern) == 1) {
        pattern        <- rep(pattern, length(c(targets$mar, targets$mnar)))
        names(pattern) <- c(targets$mar, targets$mnar)
    }
    
    ## Create a vector to hold patterns used:
    patOut        <- rep(NA, ncol(data))
    names(patOut) <- colnames(data)
    
    ## Impose MCAR missing:
    if(mechFlag["mcar"]) {
        rMat <- matrix(
            as.logical(
                rbinom(n    = prod(dim(data[ , targets$mcar])),
                       size = 1,
                       prob = pm$mcar)
            ),
            ncol = length(targets$mcar)
        )
        data[ , targets$mcar][rMat] <- NA
        patOut[targets$mcar]        <- "mcar"
    }
    
    ## Impose MAR missing:
    if(mechFlag["mar"]) {
        if(length(preds) == 1) linPred <- data[ , preds]
        else                   linPred <- rowSums(data[ , preds])
                                              
        for(i in targets$mar) {
            tmp <- makeRVec(linPred = linPred,
                            pm      = pm$mar,
                            snr     = snr$mar,
                            pattern = pattern[i])
            data[tmp$rVec, i] <- NA
            patOut[i]         <- tmp$pattern
        }
    }
    
    ## Impose MNAR missing:
    if(mechFlag["mnar"]) {
        for(i in targets$mnar) {
            tmp <- makeRVec(linPred = data[ , i],
                            pm      = pm$mnar,
                            snr     = snr$mnar,
                            pattern = pattern[i])
            data[tmp$rVec, i] <- NA
            patOut[i]         <- tmp$pattern
        }
    }
    list(data = data, pattern = patOut)
}# END imposeMissing()
