### Title:    MRE Subroutines
### Author:   Kyle M. Lang
### Created:  2019-11-13
### Modified: 2019-11-25

###--------------------------------------------------------------------------###
                                      
## Only do MI for a simple level of PM. Used for iteration planning.
mre <- function(rp, parms) {   
    ## Simulate some data:
    data <- with(parms,
                 SURF::simCovData(nObs = nObs, sigma = xCor, nVars = nVars)
                 )

    ## Run a MIBRR model:
    if(parms$doBen)
        try(
            ben(data    = data,
                y       = colnames(data)[1],
                X       = colnames(data)[-1],
                verbose = parms$verbose)
        )
    else
        try(
            bl(data    = data,
               y       = colnames(data)[1],
               X       = colnames(data)[-1],
               verbose = parms$verbose)
        )
}

###--------------------------------------------------------------------------###

## Broadcast the library function of a list of packages:
applyLib <- function(pkgList)
    lapply(pkgList, library, character.only = TRUE, logical = TRUE)

###--------------------------------------------------------------------------###

plotTrace <- function(x, y, what) {
    x <- x[[tolower(what)]]
    y <- y[[tolower(what)]]
    
    if(tolower(what) == "sigma") {
        yLim <- range(x, y)
        plot(x    = x,
             ylim = yLim,
             type = "l",
             col  = "red",
             ylab = i,
             main = paste0("Trace of ", what)
             )
        lines(y, col = "blue")
    }
    else {
        for(v in 1 : ncol(x)) {
            yLim <- range(x[ , v], y[ , v])
            plot(x    = x[ , v],
                 ylim = yLim,
                 type = "l",
                 col  = "red",
                 ylab = what,
                 main = switch(tolower(what),
                               lambda = paste0("Trace of ", what, v),
                               paste0("Trace of ",
                                      what,
                                      " for ",
                                      colnames(x)[v])
                               )
                 )
            lines(y[ , v], col = "blue")
            
            readline("Press any key to continue. ")
        }
    }
}

###--------------------------------------------------------------------------###
