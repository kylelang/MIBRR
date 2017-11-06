\name{predictMibrr}
\alias{predictMibrr}
\title{
Generate Posterior Predictions from MIBRR Models
}
\description{
  This function will generate posterior predictive draws form models fit
  using \code{miben}, \code{mibl}, \code{ben}, or \code{bl}.
}
\usage{
predictMibrr(object,
             newData,
             targetVar = NULL,
             nDraws    = 0)
}
\arguments{
  \item{object}{
    A fitted model object returned by \code{miben}, \code{mibl}, \code{ben}, or
    \code{bl}.  
  }
  \item{newData}{
    A data.frame containing new predictor data from which to generate the
    posterior predictions.   
  }
  \item{targetVar}{
    An optional character vector giving the column labels for the outcome
    variables for which to generate posterior predictions. When \code{targetVar
      = NULL} predictions are generated for all target variables contained in
    \code{object}.
  }
  \item{nDraws}{
    The number of posterior predictive draws to return.\cr
    Defaults to \code{nDraws = 0}.
  }
}
\value{
  A list containing the posterior predictive draws for each target variable
  defined in \code{targetVar} or for all target variables in \code{object} when
  \code{targetVar = NULL}.
}
\note{
  The column names of \code{newData} must contain the column names of
  all variables used to estimate the model represented by \code{object}.
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
data(predictData)

## Fit a Bayesian elastic net model:
benOut <- ben(data       = predictData$train,
              y          = "agree",
              X          = setdiff(colnames(predictData$train), "agree"),
              iterations = c(30, 10)
              )

## Generate posterior predictions for 'y':
benPred <- predictMibrr(object = benOut, newData = predictData$test)

## Fit chained Bayesian elastic net models to support MI:
mibenOut <- miben(data         = predictData$incomplete,
                  nImps        = 100,
                  iterations   = c(30, 10),
                  returnParams = TRUE)

## Generate posterior predictions from elementary imputation models:
benPred <- predictMibrr(object = mibenOut, newData = predictData$test)
}