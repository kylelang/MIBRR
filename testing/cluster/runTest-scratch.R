### Title:    Run Tests to Explore the Behavior of MCEM Chains
### Author:   Kyle M. Lang
### Created:  2019-01-01
### Modified: 2019-11-28

rm(list = ls(all = TRUE))

install.packages("devtools", repos = "http://cloud.r-project.org")
library(devtools)

install_github("kylelang/MIBRR/source/MIBRR", ref = "debugMcem")
install_github("kylelang/SURF/source/SURF", ref = "develop")

library(parallel)

verbose <- TRUE
resDir  <- "output/"
nCores  <- 2
nObs    <- 100
nVars   <- 10
pm      <- 0.0
xCor    <- 0.5
nReps   <- 8
mi      <- FALSE
sparse  <- TRUE
nPreds  <- 5
pcStart <- TRUE

source("initScript-simple.R")

si1 <- sessionInfo()
si2 <- readRDS("session_info.rds")

ls(si1)
ls(si2)

what <- ls(si2)

match <- function(x, y, what) {
    d1 <- setdiff(x[[what]], y[[what]])
    d2 <- setdiff(y[[what]], x[[what]])

    if(length(d1) + length(d2) == 0) TRUE
    else                             list("xy" = d1, "yx" = d2)
}

for(i in ls(si2)) {
    print(i)
    print(match(si1, si2, i))
}

si1
si2[7]

tmp1 <- si1[["loadedOnly"]]
tmp2 <- si2[["loadedOnly"]]
what <- names(tmp2)

check <- sapply(what, function(x) tmp1[[x]]$Version == tmp2[[x]]$Version)

what[!check]

tmp1[["survival"]]
tmp2[["survival"]]

tmp1[[what[1]]]$Version == tmp2[[what[1]]]$Version


i <- 3

what[i]
si1[[what[i]]]
si2[[what[i]]]


out <- testMcem(rp = 2, pm = pm, parms = parms, mi = mi)
                                        #class(out[[1]]$bl)
                                        #class(out[[2]]$bl)
                                        #class(out[[1]]$ben)
                                        #class(out[[2]]$ben)

## Run test in parallel:
time <- system.time(
    out <- mclapply(X        = 5 : nReps,
                    FUN      = testMcem,
                    pm       = pm,
                    parms    = parms,
                    nChains  = 2,
                    mi       = mi,
                    mc.cores = nCores)
)

out

lapply(out,
       function(y)
           lapply(y, function(x) {
               print(class(x$bl))
               print(class(x$ben))
           }
           )
       )

class(out[[1]]$bl)
class(out[[2]]$bl)
class(out[[1]]$ben)
class(out[[2]]$ben)

time / 60

saveRDS(list(mi = mi, parms = parms, out = out),
        paste0(resDir,
               "testMcemOut_",
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S"),
               ".rds")
        )
