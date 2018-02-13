\name{getImpData}
\alias{getImpData}
\title{
  Generate multiply imputed datasets.
}
\description{
  Given a fitted \code{miben} or \code{mibl} model, this funciton will
  replace missing values with imputations to generate a set of multiply
  imputed data.
}
\usage{
getImpData(mibrrFit, nImps)
}
\arguments{
  \item{mibrrFit}{
    A fitted model object (with class \code{MibrrFit}) returned by
    \code{miben}, \code{mibl}.  
  }
  \item{nImps}{
    The number of imputed datasets to return.
  }
}
\value{
  A list containing \code{nImps} imputed datasets wherein the missing
  values in the target variables named in the \code{mibrrFit} object are
  replaced by random draws from the appropriate posterior predicted
  distributions.
}
\note{
  \code{nImps} cannot be larger than the posterior Gibbs sample size in the
  \code{mibrrFit} object.
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
\code{\link{miben}}, \code{\link{mibl}}
}
\examples{
data(mibrrExampleData)

## Fit a Bayesian elastic net model:
mibenOut <- miben(data         = mibrrExampleData,
                  targetVars   = c("y", paste0("x", c(1 : 3))),
                  ignoreVars   = "idNum",
                  sampleSizes  = c(500, 500),
                  doMcem       = FALSE,
                  lam1PriorPar = c(1.0, 0.1),
                  lam2PriorPar = c(1.0, 0.1)
                  )

## Generate 25 imputed datasets:
mibenImps <- getImpData(mibrrFit = mibenOut, nImps = 25)
}