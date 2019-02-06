\name{ppCheck}
\alias{ppCheck}
\title{
Posterior Predictive Checks
}
\description{
  This function will generate overlaid density plots of the posterior
  predictive distributions and the observed outcomes.
}
\usage{
ppCheck(mibrrFit, targetVars = NULL)
}
\arguments{
  \item{mibrrFit}{
    An estimated \code{MibrrFit} object.
  }
  \item{targetVars}{
    A character vector naming the variables to check. When
    \code{targetVars = NULL}, all variables listed as \code{targetVars}
    in the object provided to \code{mibrrFit} will be checked.  
  }
}
\value{
  Used for its side-effects.
}
\author{
  Kyle M. Lang
}
\section{Warning}{
  This function is in a highly unstable \emph{alpha} level of development. Please
  anticipate frequent---and dramatic---changes to the functionality and user
  interface.
  
  You have been granted access to this package for evaluation purposes, only. This
  function is \strong{absolutely not} ready for use in real-world analyses!
}
\seealso{
\code{\link{getParams}}, \code{\link{getField}}
}
\examples{
data(mibrrExampleData)

## Estimate a Bayesian linear regression model:
bvrOut <- bvr(data        = mibrrExampleData,
              y           = "y",
              X           = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
              sampleSizes = c(500, 500),
              control     = list(savePpSams = TRUE)
              )

## Generate PP Plots for "y":
ppCheck(bvrOut)
}