### Title:    Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2015-01-01
### Modified: 2019-11-15

rm(list = ls(all = TRUE))

outDir    <- "../output/"

files <- system(paste0("ls ", outDir), intern = TRUE)
files <- files[grep("^testMcemOut", files)]

                                        #tmp <- readRDS(paste0(outDir, files[length(files)]))

tmp <- readRDS(paste0(outDir, files[5]))

out   <- tmp$out
parms <- tmp$parms

parms

tmp <- out[[i]][[1]]$bl
ls(tmp)

tmp$iterations

## MIBL:
i <- 7
par(mfcol = c(1, 1))
l1 <- out[[i]][[1]]$bl$lambdaHistory[[1]][ , 1]
l2 <- out[[i]][[2]]$bl$lambdaHistory[[1]][ , 1]

plot(x    = l1,
     ylim = range(l1, l2),
     type = "l",
     main = "Lambda",
     col  = "red")
lines(l2, col = "blue")


n   <- parms$iters["nEmTune"] / 2
l12 <- l1[(length(l1) - n + 1) : length(l1)]
l22 <- l2[(length(l2) - n + 1) : length(l2)]

tmp <- matrix(c(l12, l22), ncol = 2)

b <- mean(apply(tmp, 1, var))
w <- mean(apply(tmp, 2, var))

b / w

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
