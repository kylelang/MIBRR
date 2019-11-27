### Title:    Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2015-01-01
### Modified: 2019-11-15

rm(list = ls(all = TRUE))

outDir    <- "../output/"

files <- system(paste0("ls ", outDir), intern = TRUE)
files <- files[grep("^testMcemOut", files)]

tmp <- readRDS(paste0(outDir, files[length(files)]))

out   <- tmp$out
parms <- tmp$parms

parms

## Find failed reps:
check <- c()
for(i in 1 : length(out))
    check[i] <- class(out[[i]][[1]]$ben) == "try-error" |
        class(out[[i]][[2]]$ben) == "try-error"

which(check)
goodReps <- setdiff(1 : length(out), which(check))
goodReps

### Print Stuff ###

## MIBEN:
par(mfcol = c(1, 2))
for(i in goodReps) {
    l1 <- out[[i]][[1]]$ben$lambdaHistory[[1]][ , 1]
    l2 <- out[[i]][[2]]$ben$lambdaHistory[[1]][ , 1]
    
    plot(x    = l1,
         ylim = range(l1, l2),
         type = "l",
         main = "L1",
         col  = "red")
    lines(l2, col = "blue")
    
    l1 <- out[[i]][[1]]$ben$lambdaHistory[[1]][ , 2]
    l2 <- out[[i]][[2]]$ben$lambdaHistory[[1]][ , 2]
    
    plot(x    = l1,
         ylim = range(l1, l2),
         type = "l",
         main = "L2",
         col  = "red")
    lines(l2, col = "blue")
    
    readline("Hit any key to continue. ")
}

## Find failed reps:
check <- c()
for(i in 1 : length(out))
    check[i] <- class(out[[i]][[1]]$bl) == "try-error" |
        class(out[[i]][[2]]$bl) == "try-error"

which(check)
goodReps <- setdiff(1 : length(out), which(check))
goodReps

## MIBL:
par(mfcol = c(1, 1))
for(i in goodReps) {
    l1 <- out[[i]][[1]]$bl$lambdaHistory[[1]][ , 1]
    l2 <- out[[i]][[2]]$bl$lambdaHistory[[1]][ , 1]
    
    plot(x    = l1,
         ylim = range(l1, l2),
         type = "l",
         main = "Lambda",
         col  = "red")
    lines(l2, col = "blue")

    readline("Hit any key to continue. ")
}

## Check error messages:
for(i in which(check)) {
    tmp <- out[[i]][[1]]$bl
    print(tmp)
}

## Markov Chains ##

## MIBEN:
par(mfrow = c(1, 1))
for(i in goodReps) {
    tmp1 <- getParams(out[[i]][[1]]$ben, 1)
    tmp2 <- getParams(out[[i]][[2]]$ben, 1)
    
    plotTrace(tmp1, tmp2, "Beta")
    plotTrace(tmp1, tmp2, "Tau")
    plotTrace(tmp1, tmp2, "Sigma")
    
    readline("Hit any key to continue. ")
}

## MIBL:
for(i in 1 : length(out)) {
    tmp1 <- getParams(out[[i]][[1]]$bl, 1)
    tmp2 <- getParams(out[[i]][[2]]$bl, 1)
    
    plotTrace(tmp1, tmp2, "Beta")
    plotTrace(tmp1, tmp2, "Tau")
    plotTrace(tmp1, tmp2, "Sigma")

    readline("Hit any key to continue. ")
}
