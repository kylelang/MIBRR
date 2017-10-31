\name{miben}
\alias{miben}
\title{
Multiple imputation with the Bayesian elastic net.
}
\description{
  This function implements MIBEN, a robust MICE-based multiple imputation scheme
  that employs the Bayesian elastic net as its elementary imputation method.
}
\usage{
miben(data,
      nImps          = 100,
      targetVars     = NULL,
      ignoreVars     = NULL,
      iterations     = c(100, 10),
      sampleSizes    = list(rep(25, 2), rep(250, 2), rep(500, 2)),
      missCode       = NULL,
      returnConvInfo = TRUE,
      returnParams   = FALSE,
      verbose        = TRUE,
      seed           = NULL,
      control)  
}
\arguments{
  \item{data}{
    An incomplete, numeric data matrix or data frame for which to create the
    imputations.
  }
  \item{nImps}{
    An integer giving the number of imputations to create. Defaults to
    \code{nImps = 100}.
  }
  \item{targetVars}{
    An optional character vector giving the column labels for the variables to
    be imputed. When \code{targetVars = NULL}, all variables not listed in
    \code{ignoreVars} are imputed.  
  }
  \item{ignoreVars}{
    An optional character vector giving the column labels for those variables
    that should be excluded from the imputation model (e.g., ID variables).
  }
  \item{iterations}{
    A two-element numeric vector giving the number of iterations to employ
    during the MCEM approximation and tuning phases, respectively. Defaults to
    \code{iterations = c(100, 10)}.
  }
  \item{sampleSizes}{
    A list containing three two-element numeric vectors giving the number of
    MCMC draws discarded as burn-in and retained, respectively, during the MCEM
    approximation, tuning, and sampling phases. Defaults to \code{sampleSizes =
      list(rep(25, 2), rep(250, 2), rep(500, 2))}. 
  }
  \item{missCode}{
    An optional integer-valued code used to flag the missing data in
    \code{rawData}. Should take a value that cannot naturally occur in
    \code{rawData}. Not needed when the missing data are coded as \code{NA}.
  }
  \item{returnConvInfo}{
    A logical switch: should convergence information for the imputation model
    (i.e., history of the optimized penalty parameters and R-Hat values for
    final parameter estimates) be returned? Defaults to \code{returnConvInfo =
      TRUE}. 
  }
  \item{returnParams}{
    A logical switch: should the final Gibbs samples of the imputation model's
    parameters be returned? Defaults to \code{returnParams = FALSE}.
  }
  \item{verbose}{
    A logical switch: should verbose output be printed to stdout? Defaults to
    \code{verbose = TRUE}.
  }
  \item{seed}{
    An integer-valued seed for the pseudorandom number generator. When
    \code{seed = NULL} R's default PRNG and seed are left alone.
  }
  \item{control}{
    A list of control parameters for the Gibbs sampler and penalty parameter
    optimization (see Details for more information).
  }
}
\details{
  \code{control} is a list containing six elements:
  \itemize{
    \item convThresh = The R-Hat value used to judge convergence. R-Hat
    values < convThresh arising during the MCEM tuning phase will trigger a
    warning. Defaults to \code{convThres = 1.1}. 
    \item lambda1Starts = An optional numeric vector giving starting values for
    the LASSO penalty parameter. Values are recycled to populate a vector with
    \code{size = length(targetVars)}. Defaults to \code{rep(0.5,
      length(targetVars))}.
    \item lambda2Starts = An optional numeric vector giving starting values for
    the ridge penalty parameter. Values are recycled to populate a vector with
    \code{size = length(targetVars)}. Defaults to \code{rep(0.1 * nPreds,
      length(targetVars))}, where nPreds is the number of predictors in the
    imputation model.
    \item usePcStarts = A logical switch: use the starting values for lambda1
    suggested by Park and Casella (2009)? Defaults to \code{usePcStarts = FALSE}.
    \item smoothingWindow = An integer giving the number of approximation phase
    lambda values to average over to get the starting lambdas for the MCEM
    tuning phase. Setting smoothingWindow > 1 can facilitate convergence of the
    MCEM tuning phase when burn-in lambda estimates are very noisy. Defaults to
    \code{smoothingWindow = 1}.
    \item center = A logical switch: Should the data be centered before
    estimating the imputation model? When \code{center = TRUE} imputed data are
    un-centered when returned. Defaults to \code{center = TRUE}.
    \item scale = A logical switch: Should the predictor data be scaled to have
    unit variance before estimating the imputation model? When \code{scale =
      TRUE} imputed data are un-scaled when returned. Defaults to \code{scale =
      TRUE}.
    \item adaptScales = A logical switch: Should the target variables' scales be
    actively updated as part of imputation model estimation? Defaults to
    \code{adaptScales = TRUE}.
    \item simpleIntercept
    \item minPredCor
    \item miceIters
    \item miceMethod
    \item fimlStarts
    \item preserveStructure
    \item optTraceLevel
    \item optCheckKkt
    \item optMethod
    \item optBoundLambda
    \item optReturnAcov
  }
}
\value{
  A list containing \code{nImps} imputed versions of \code{data} and whatever
  additional output is requested via \code{returnConvInfo} and
  \code{returnParams}.
}
\references{
  Lang, K. M. (2015). \emph{MIBEN: Multiple imputation with the
    Bayesian elastic net} (Unpublished doctoral dissertation). University of
  Kansas.
  
  Li, Q., & Lin, N. (2010). The Bayesian Elastic Net. \emph{Bayesian
    Analysis, 5}(1), 151--170.
}
\author{
  Kyle M. Lang

  Maintainer: Kyle M. Lang <k.m.lang@uvt.nl>
}
\note{
You are not a wombat...or are you???
}
\seealso{
\code{\link{mibl}}
}
\examples{
data(mibrrExampleData)

miben(data       = mibrrExampleData,
      targetVars = c("y", paste0("x", c(1 : 3))),
      ignoreVars = "idNum")
}