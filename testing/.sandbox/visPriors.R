### Title:    Visualize Penalty Parameters' Priors
### Author:   Kyle M. Lang
### Created:  2019-MAR-06
### Modified: 2019-MAR-06

### BL ###
shape <- 10.1
rate  <- 0.1

mn <- shape / rate
vr <- shape / rate^2

mn
vr

end    <- 1000
badSam <- TRUE
while(badSam) {
    x <- seq(0, end, 0.1)
    y <- dgamma(x, shape = shape, rate = rate)
    
    ## Filter out trivial density values:
    filter <- y > 1e-5

    ## Do we have adequate support?
    if(filter[end])
        end <- end + 100
    else 
        badSam <- FALSE
}

## Hack to keep initial element in highly skewed distributions:
tmp <- min(which(filter))
filter[tmp - 1] <- TRUE

plot(x = x[filter], y = y[filter], type = "l")
