\name{inspect}
\alias{inspect}
\title{
  Access arbitrary fields in a MibrrFit object.
}
\description{
  This function will access any field in a fitted \code{MibrrFit} object
  using a familiar, S3-like interface.
}
\usage{
inspect(mibrrFit, what)
}
\arguments{
  \item{mibrrFit}{
    A fitted model object (with class \code{MibrrFit}) returned by
    \code{miben}, \code{mibl}.  
  }
  \item{what}{
    The name of the field to extract (see Details for a list of the most
    interesting possibilities).
  }
}
\value{
  The contents of whatever field is defined by the \code{what}
  argument.
}
\details{
  This function can access any field in the \code{MibrrFit} object (see
  \code{\link{MibrrFit}} for a full list of possible fields), but the
  most interesting possibilities are the following:

  \describe{
    \item{data}{}
    \item{rHats}{}
    \item{lambdaHistory}{}
    \item{gibbsOut}{}
  }
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
              ignoreVars   = "idNum",
              sampleSizes  = c(500, 500),
              doMcem       = FALSE,
              lam1PriorPar = c(1.0, 0.1)
              lam2PriorPar = c(1.0, 0.1)
              )

## Extract the potential scale reduction factors:
mibenPars <- inspect(mibrrFit = benOut, what = "rHats")
}