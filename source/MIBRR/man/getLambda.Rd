\name{getLambda}
\alias{getLambda}
\title{
  Extract samples of penalty parameters.
}
\description{
  Given a fitted \code{miben}, \code{mibl}, \code{ben}, or \code{bl}
  model, this funciton will extract the stationary samples of the
  penalty parameters.
}
\usage{
getLambda(mibrrFit, target)
}
\arguments{
  \item{mibrrFit}{
    A fitted model object (with class \code{MibrrFit}) returned by
    \code{miben}, \code{mibl}.  
  }
  \item{target}{
    The column label for the target variable whose parameter samples are
    to be extracted.
  }
}
\value{
  A list containing the stationary samples of the penalty parameters
  for the model wherein \code{target} was the outcome variable.
}
\note{
  If \code{doMcem = TRUE} the final optimized values of the penalty
  parameters, \eqn{\{\lambda_1, \lambda_2\}}{\{lambda1, lambda2\}}, are
  returned. If \code{doMcem = FALSE}, the posterior samples of the
  penalty parameters are returned.
}
\references{
  Lang, K. M. (2015) \emph{MIBEN: Multiple imputation with the
    Bayesian elastic net} (Unpublished doctoral dissertation). University of
  Kansas.
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
\code{\link{miben}}, \code{\link{mibl}}, \code{\link{ben}}, \code{\link{bl}}
}
\examples{
data(mibrrExampleData)

## Fit a Bayesian elastic net model:
benOut <- ben(data         = mibrrExampleData,
              y            = "y",
              X            = paste0("x", c(1 : 3)),
              sampleSizes  = c(500, 500),
              doMcem       = FALSE,
              lam1PriorPar = c(1.0, 0.1),
              lam2PriorPar = c(1.0, 0.1)
              )

## Extract penalty parameter samples:
lambdas <- getLambda(mibrrFit = benOut, target = "y")
}