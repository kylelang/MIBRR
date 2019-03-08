### Title:    Visualize Penalty Parameters' Priors
### Author:   Kyle M. Lang
### Created:  2019-MAR-06
### Modified: 2019-MAR-06

library(HyperbolicDist)

visPriors <- function(p1, p2, what, eps = 1e-4) {
    ## Define the range of x:
    supp <- switch(what,
                   l1 = qgamma(c(eps, 1 - eps), shape = p1, rate = p2),
                   l2 = qgig(c(eps, 1 - eps), Theta = c(1.0, p1, p2))
                   )

    ## Generate the density points:
    x <- seq(supp[1], supp[2], length.out = 1000)
    y <- switch(what,
                l1 = dgamma(x, shape = p1, rate = p2),
                l2 = dgig(x, Theta = c(1.0, p1, p2))
                )
    
    ## Plot the prior density:
    plot(x = x,
         y = y,,
         type = "l",
         ylab = "Density",
         xlab = switch(what, l1 = "Lambda1", l2 = "Lambda2^2"),
         main = paste0("Prior for ",
                       switch(what, l1 = "Lambda1", l2 = "Lambda2^2"),
                       "\n",
                       switch(what, l1 = "Shape = ", l2 = "Chi = "),
                       p1,
                       "\n",
                       switch(what, l1 = "Rate = ", l2 = "Psi = "),
                       p2)
         )
}

par(mfrow = c(5, 5))

s1 <- seq(1, 100, length.out = 5)
s2 <- seq(0.1, 0.001, length.out = 5)

for(i in s1)
    for(j in s2)
        visPriors(i, j, "l2")
v <- 0.1
par(mfrow = c(4, 4))
for(m in seq(0, 2, length.out = 16)[-1])
    visPriors(m^2/v, m/v, "l1")

visPriors(1, 1, "l1")
visPriors(0.5, 5.0, "l2")

(1 : 16)^2 / (1 : 16)

1 / 1.5
1 / 1.78
sqrt(.237)

library(MIBRR)

?miben
