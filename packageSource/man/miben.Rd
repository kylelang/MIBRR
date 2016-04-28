\name{miben}
\alias{miben}
\title{
Multiple imputation with the Bayesian elastic net.
}
\description{
  This function implements MIBEN, a robust MICE-based multiple imputation
  scheme that employs the Bayesian elastic net as its elementary imputation
  method.
}
\usage{
miben(rawData, nImps = 100, targetVariables = NULL,
      ignoreVariables = NULL, nEmApproxIters = 100, nEmTuneIters = 10,
      nEmApproxBurn = 25, nEmApproxGibbs = 25, nEmTuneBurn = 100,
      nEmTuneGibbs = 200, nPosteriorBurn = NULL, nPosteriorThin = 5,
      missCode = NULL, returnConvInfo = TRUE, returnModelParams = FALSE,
      verboseIters = TRUE, verboseErrors = TRUE, gibbsControl,
      optControl)  
}
\arguments{
  \item{rawData}{
    an incomplete, numeric data matrix or data frame for which to create
    the imputations
  }
  \item{nImps}{
    the number of imputations to create. Defaults to 100.
  }
  \item{targetVariables}{
    an optional character vector giving the column labels for the
    variables to be imputed. Defaults to imputing all variables not
    listed in \code{ignoreVariables} if left blank.  
  }
  \item{ignoreVariables}{
    an optional character vector giving the column labels for those variables
    that should be excluded from the imputation model (e.g., ID variables).
  }
  \item{nEmApprox}{
    an integer giving the number of iterations to employ during the
    initial MCEM approximation phase. Defaults to 100. 
  }
  \item{nEmTune}{
    an integer giving the number of iterations to employ during the MCEM
    tuning phase. Defaults to 10. 
  }
  \item{emApproxBurn}{
    an integer giving the number of "burn in" Markov chain iterations to
    be discarded during each Gibbs sampling iteration of the MCEM
    approximiation phase. Defaults to 25.
  }
  \item{emApproxGibbs}{
    an integer giving the size of the retained Gibbs samples drawn during
    the MCEM approximation phase. Defaults to 25.
  }
  \item{emTuneBurn}{
    an integer giving the number of "burn in" Markov chain iterations to
    be discarded during each Gibbs sampling iteration of the MCEM
    tuning phase. Defaults to 100. 
  }
  \item{emTuneGibbs}{
    an integer giving the size of the retained Gibbs samples drawn
    during the tuning phase of the MCEM algorithm. Defaults to 200.
  }
  \item{posteriorBurn}{
    an optional integer giving the number of burn-in Markov chain
    iterations to discard when constructing the stationary posterior
    distribution. Will default to \code{emTuneBurn} if left empty.
  }
  \item{posteriorThin}{
    an integer giving the thinning factor used to construct the final
    posterior distribution. \code{posteriorThin = 1} indicates no
    thinning of the posterior predictive draws. Defaults to 5 (see
    'Details' for more information).
  }
  \item{missCode}{
    an optional integer-valued code used to flag the missing data in
    \code{rawData}. Should take a value that cannot naturally occur
    in \code{rawData}. Not needed when the missing data are coded as
    \code{NA}.
  }
  \item{returnConvInfo}{
    a logical switch, should convergence information for the imputation
    model (i.e., history of the optimized penalty parameters and R-Hat
    values for final parameter estimates) be returned? Defaults to
    TRUE.
  }
  \item{returnModelParams}{
    a logical switch, should the final Gibbs samples of the imputation
    model's parameters be returned? Defaults to FALSE.
  }
  \item{verboseIters}{
    a logical switch, should iteration history be printed to stdout?
    Defaults to TRUE.
  }
  \item{verboseErrors}{
    a logical switch, in the event of an error, should additional
    diagnostic information be printed to stderr? Defaults to TRUE.
  }
  \item{gibbsControl}{
    a list of control parameters for the Gibbs sampler (see Details for
    more information)
  }
  \item{optControl}{
    a list of control parameters for the optimization of the penalty
    parameters (see Details for more information)
  }
}
\details{
  \code{posteriorThin} dictates the number of draws used to approximate
  the posterior predictive distribution of the missing data data
  according to the rule: \code{size = posteriorThin * nImps}. The actual
  imputations are produced by randomly sampling (without replacement)
  from this approximated posterior predictive distribution. Thus,
  \code{posteriorThin} only induces \emph{thinning} in an indirect,
  asymptotic sense.
  
  \code{gibbsControl} is a list containing six elements:
  \itemize{
    \item createRngStream = a logical switch indicating whether a new
    L'ecuyer RNG stream should be created. Must equal TRUE unless an
    extant stream name is given as \code{streamName}. Defaults to TRUE. 
    \item streamName = an optional character string giving the stream
    name corresponding to an extant L'ecuyer RNG stream such as one
    created by the package \code{rlecuyer}.
    \item convThresh = the R-Hat value used to judge
    convergence. R-Hat values < convThresh arising during the MCEM
    tuning phase will return a warning. Defaults to 1.1. 
    \item rngSeed = an integer used to seed the L'ecuyer RNG
    stream. Defaults to 235711.
  }
  
  \code{optControl} is a list containing four elements:
  \itemize{
    \item lambda1Starts = an optional numeric vector giving starting
    values for the LASSO penalty parameter. Values are recycled to
    populate a vector with \code{size = length(targetVars)}. Defaults to
    \code{rep(0.5, length(targetVars))}, if left empty.
    \item lambda2Starts = an optional numeric vector giving starting
    values for the ridge penalty parameter. Values are recycled to
    populate a vector with \code{size = length(targetVars)}. Defaults to
    \code{rep(0.1 * nPreds, length(targetVars))}, if left empty where
    nPreds is the number of predictors in the imputation model.
    \item emEpsilon = a numeric value giving the convergence criterion
    used to judge convergence of the marginal, numeric optimization of
    the penalty parameters. This value is acts as an nlopt stopping rule
    and indicates the threshold for the relative change in the function
    value corresponding to convergence. Defaults to 1.0e-5. 
    \item smoothingWindow = an integer giving the number of burn-in
    lambda values to average over to get the starting lambdas for the
    MCEM tuning phase. Setting smoothingWindow > 1 can facilitate
    convergence of the MCEM tuning phase when burn-in lambda estimates
    are very noisy. Defaults to 1.
  }
}
\value{
  A list containing \code{nImps} imputed versions of \code{rawData} and
  whatever additional output is requested via \code{returnConvInfo} and
  \code{returnModelParams}. 
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

  Maintainer: Kyle M. Lang <kylelang@ku.edu>
}
\note{
You are not a wombat...or are you???
}
\seealso{
\code{\link{mibl}}
}
\examples{
data(mibrrExampleData)

miben(rawData = mibrrExampleData,
      targetVars = c("y", paste0("x", c(1 : 3))),
      ignoreVars = "idNum")
}