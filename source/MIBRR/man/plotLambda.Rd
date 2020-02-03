\name{plotLambda}
\alias{plotLambda}
\title{
  Plot lambda's MCEM chains
}
\description{
  Given a fitted \code{miben}, \code{mibl}, \code{ben}, or \code{bl}
  model, this funciton will plot (overlaid) tracelines of each penalty
  parameter's (or the loglikelihood's) MCEM chains.
}
\usage{
plotLambda(mibrrFit, target, logLik = FALSE)
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
  \item{logLik}{
    A logical switch: Plot the model's loglikelihood value instead of
    the penalty parameter(s) value(s)?
  }
}
\value{
  None. This function is used for its side-effects.
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

## Fit a Bayesian LASSO model:
blOut <- bl(data    = mibrrExampleData,
            y       = "y",
            X       = paste0("x", c(1 : 3)),
            nChains = 2)

## Plot lambda samples:
plotLambda(mibrrFit = blOut, target = "y")
}