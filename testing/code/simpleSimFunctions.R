

cleanData <- function(parms) {
    dat1 <- readRDS(paste0(parms$dataDir, "adamsKlpsData.rds"))
    
    sysRacNames <- paste0("NORI", c(2, 3, 6, 7, 8, 9, 11))
    sysRacNames <- sysRacNames[sysRacNames %in% colnames(dat1)]
    
    indRacNames <- paste0("NORI", c(1, 4, 5, 10))
    indRacNames <- indRacNames[indRacNames %in% colnames(dat1)]
    
    policyNames <- paste0("POLICY", c(1, 3, 4, 5, 6))
    policyNames <- policyNames[policyNames %in% colnames(dat1)]
    
    meritNames <- paste0("WPRIV", c(1 : 10))
    meritNames <- meritNames[meritNames %in% colnames(dat1)]

    data.frame(
        sysRac = rowMeans(dat1[ , sysRacNames], na.rm = TRUE),
        indRacData = rowMeans(dat1[ , indRacNames], na.rm = TRUE),
        policyData = rowMeans(dat1[ , policyNames], na.rm = TRUE),
        meritData = rowMeans(dat1[ , meritNames], na.rm = TRUE)
    )
}


cleanData2 <- function(parms) {
    dat1 <- readRDS(paste0(parms$dataDir, "adamsKlpsData.rds"))
    dat1[ , c("sysRac", "indRac", "policy", "merit", "polAffil", "revDisc")]
}


simSimpleData <- function(nObs, r2, alpha, beta, collinVal) {
    sigma <- matrix(collinVal, 3, 3)
    diag(sigma) <- 1.0
    
    X <- rmvnorm(nObs, parms$predMeans, sigma)
    eta <- X %*% beta
    errorVar <- (var(eta) / r2) - var(eta)
    y <- alpha + eta + rnorm(nObs, 0, sqrt(errorVar))
    outData <- data.frame(y, X)
    colnames(outData) <- c("y", "x1", "x2", "z")
    outData
}


simData2 <- function(nObs, r2, alpha, beta, collinVal) {
    sigma <- matrix(collinVal, 3, 3)
    diag(sigma) <- 1.0
    
    X <- rmvnorm(nObs, rep(0, 3), sigma)
    eta <- X[ , -3] %*% beta
    errorVar <- (var(eta) / r2) - var(eta)
    y <- alpha + eta + rnorm(nObs, 0, sqrt(errorVar))
    outData <- data.frame(y, X)
    colnames(outData) <- c("y", "x1", "x2", "z")
    outData
}


simData3 <- function(parms) {
    as.data.frame(rmvnorm(parms$nObs, parms$mu, parms$sigma))
}


imposeMissing <- function(inData, parms)
{
    incompVars <- parms$incompVars
    auxVars <- parms$auxVars
    pm <- parms$pm
    marType <- vector("character", length(incompVars))
    marType <- parms$marType
    
    if(length(auxVars) > 1) auxVar <- rowMeans(inData[ , auxVars])
    else                    auxVar <- inData[ , auxVars]
    
    pVec <- pnorm(auxVar, mean(auxVar), sd(auxVar))
    
    for(v in 1 : length(incompVars)) {
        if(marType[v] == "tails") {
            rVec <- pVec < (pm / 2) | pVec > (1 - (pm/2))
        } else if(marType[v] == "center") {
            rVec <- pVec > (0.5 - (pm / 2)) & pVec < (0.5 + (pm / 2))
        } else if(marType[v] == "lower") {
            rVec <- pVec < pm
        } else if(marType[v] == "upper") {
            rVec <- pVec > (1 - pm)
        } else {
            stop("Please provide a valid 'marType'")
        }
        inData[rVec, incompVars[v]] <- NA
    }
    inData
}


doRep <- function(rp, parms) {
    cat(paste0("Doing replication ", rp, "\n"))
    testForm <- parms$testForm
    incompVars <- parms$incompVars
    resDir <- parms$resDir
    resList <- list()
    
    dat1 <- simData3(parms = parms)
    
    resList$comp <- lm(testForm, data = dat1)
    
    missData <- imposeMissing(inData = dat1, parms = parms)
    
    suppressWarnings(
        mibenOut <- miben(rawData         = missData,
                          targetVariables = incompVars,
                          verboseIters    = FALSE,
                          verboseErrors   = FALSE)
    )
    
    fitList <- lapply(mibenOut$imps, FUN = function(x) lm(testForm, data = x))
    resList$miben <- MIcombine(fitList)
    
    miceOut <- mice(data      = missData,
                    m         = 100,
                    maxit     = 5,
                    printFlag = FALSE)
    
    impList <- list()
    for(m in 1 : 100) impList[[m]] <- complete(miceOut, m)
    
    fitList <- lapply(impList, FUN = function(x) lm(testForm, data = x))
    resList$mice <- MIcombine(fitList)

    saveRDS(resList, paste0(resDir, "outFile", rp, ".rds"))
    list()# Return nothing
}
