### Title:    Testing Subroutines
### Author:   Kyle M. Lang
### Created:  2015-JAN-01
### Modified: 2019-FEB-26

###--------------------------------------------------------------------------###

mcem0 <- function(y,
                  X,
                  l2,
                  s2,
                  beta,
                  nIter,
                  sams,
                  intercept = TRUE,
                  norm      = TRUE)
{
    ## Initialize lambda's history:
    l2 <- c(l2, rep(NA, (nIter - 1)))
    
    for(i in 1 : nIter) {
        ## Updated Gibbs samples:
        out0 <- blasso(y         = y,
                       X         = X,
                       T         = sams[2],
                       thin      = sams[1],
                       lambda2   = l2[i],
                       s2        = s2,
                       beta      = beta,
                       rd        = FALSE,
                       RJ        = FALSE,
                       rao.s2    = FALSE,
                       icept     = intercept,
                       normalize = norm)
        
        ## Update parameter starting values:
        s2   <- calcMode(out0$s2, discrete = FALSE)
        beta <- apply(out0$beta, 2, calcMode, discrete = FALSE)
        
        ## Optimize lambda:
        if(i < nIter)
            l2[i + 1] <-
                2 * ncol(X) / sum(colMeans(1 / out0$tau2i, na.rm = TRUE))
    }
    
    list(lambda = l2, s2 = s2, beta = beta)
}

###--------------------------------------------------------------------------###

bl0Mcem <- function(data,
                    yName,
                    xNames,
                    iters,
                    sams,
                    l2        = 1.0,
                    intercept = TRUE,
                    norm      = TRUE)
{
    ## Partition data:
    y <- data[ , yName]
    X <- data[ , xNames]
    
    ## Initialize parameter starting values:
    s2     <- var(y - mean(y))
    beta   <- rnorm(ncol(X))
    lambda <- c()
    
    ## Run approximation and tuning iterations:
    for(s in 1 : 2) {
        ## Run MCEM with a given specification of {nIter, sams}:
        out <- mcem0(y         = y,
                     X         = X,
                     l2        = l2,
                     s2        = s2,
                     beta      = beta,
                     nIter     = iters[s],
                     sams      = sams[[s]],
                     intercept = intercept,
                     norm      = norm)
        
        ## Update parameter starting values:
        l2   <- out$lambda[iters[s]]
        s2   <- out$s2
        beta <- out$beta
        
        ## Extend lambda's history:
        lambda <- c(lambda, out$lambda)
    }

    ## Run final model with optimized lambda:
    out <- blasso(y       = y,
                  X       = X,
                  T       = sams[[3]][2],
                  thin    = sams[[3]][1],
                  lambda2 = l2,
                  s2      = s2,
                  beta    = beta,
                  rd      = FALSE,
                  RJ      = FALSE,
                  rao.s2  = FALSE,
                  icept   = intercept)

    list(out = out, lambda = lambda)
}

###--------------------------------------------------------------------------###

predBl0 <- function(obj, X, type = "raw") {
    xb <- X %*% t(obj$beta)
    e  <- rnorm(n = length(obj$s2), sd = sqrt(obj$s2))
    
    out <- t(t(xb) + obj$mu + e)
    
    if(type == "eap") out <- rowMeans(out) 
    if(type == "map") out <- apply(out, 1, calcMode, discrete = FALSE)
    
    out
}

###--------------------------------------------------------------------------###

predBen0 <- function(obj, X, type = "raw") {
    ## Extract parameter samples:
    tmp <- extract(obj, pars = c("mu", "beta", "sigma2"))
    
    xb <- X %*% t(tmp[["beta"]])
    b0 <- as.numeric(tmp[["mu"]])
    
    s2 <- tmp[["sigma2"]]
    e  <- rnorm(n = length(s2), sd = sqrt(s2))
    
    out <- t(t(xb) + b0 + e)
    
    if(type == "eap") out <- rowMeans(out) 
    if(type == "map") out <- apply(out, 1, calcMode, discrete = FALSE)
    
    out
}

###--------------------------------------------------------------------------###

madOutliers <- function(x, cut = 2.5, na.rm = TRUE) {
    ## Compute the median and MAD of x:
    mX   <- median(x, na.rm = na.rm)
    madX <- mad(x, na.rm = na.rm)
    
    ## Return row indices of observations for which |T_MAD| > cut:
    which(abs(x - mX) / madX > cut)
} 

###--------------------------------------------------------------------------###

## Define a simple Bayesian regression function:
bReg <- function(data, y, X, nSams, scale = "none") {
    ## Define the model formula:
    f1    <- paste(y, paste(X, collapse = " + "), sep = " ~ ")
    
    if(scale == "X")
        data[ , X] <- scale(data[ , X])
    else if(scale == "y")
        data[ , y] <- scale(data[ , y])
    else if(scale == "all")
        data <- scale(data)
    
    ## Get the expected betas via least squares:
    fit   <- lm(f1, data = data)
    beta0 <- coef(fit)
    
    ## Sample sigma:
    sigma2 <- rinvchisq(nSams, df = fit$df, scale = summary(fit)$sigma^2)
    
    ## Sample beta:
    beta           <- matrix(NA, nSams, length(X) + 1)
    colnames(beta) <- paste0("b", 0 : length(X))
    for(n in 1 : nSams) {
        betaVar   <- sigma2[n] * solve(crossprod(qr.X(fit$qr)))
        beta[n, ] <- rmvnorm(1, mean = beta0, sigma = betaVar)
    }
    
    ## Return the posterior samples:
    list(beta = beta, sigma2 = sigma2)
}
