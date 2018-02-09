\name{ben}
\alias{ben}
\title{
Bayesian Elastic Net
}
\description{
  This function will fit the Bayesian elastic net to (possibly) incomplete data.
}
\usage{
ben(data,
    y,
    X            = NULL,
    iterations   = c(100, 10),
    sampleSizes  = list(rep(25, 2), rep(250, 2), rep(500, 2)),
    doMcem       = TRUE,
    lam1PriorPar = NULL,
    lam2PriorPar = NULL,
    missCode     = NULL,
    verbose      = TRUE,
    seed         = NULL,
    control      = list()
    )  
}
\arguments{
  \item{data}{
    A, possibly incomplete, numeric data matrix or data frame to which the BEN
    is to be fit.
  }
  \item{y}{
    The column label for the outcome variable.  
  }
  \item{X}{
    An optional character vector giving the column labels for the predictor
    variables. When \code{X = NULL} the target variable is regressed onto all
    other variables in \code{data}.
  }
  \item{iterations}{
    A two-element numeric vector giving the number of iterations to employ
    during the MCEM approximation and tuning phases, respectively.\cr
    Defaults to \code{iterations = c(100, 10)}.
  }
  \item{sampleSizes}{
    A list or vector giving the desired Gibbs sample sizes (see
    Details for more information).\cr 
    Defaults to \code{sampleSizes = list(rep(25, 2), rep(250, 2),
      rep(500, 2))}. 
  }
  \item{doMcem}{
    A logical switch: Should the model be estimated using Markov Chain
    Expectation Maximization (MCEM)? If \code{doMcem = FALSE} the model
    is estimated as a fully Bayesian model by assigning hyper-priors to
    the penalty parameters.\cr
    Defaults to \code{doMcem = TRUE}.
  }
  \item{lam1PriorPar}{
    A two-element numeric vector giving the prior shape and rate
    parameters, respectively, for the squared LASSO penalty parameter's,
    \eqn{\lambda_1^2}{lambda1^2}, gamma prior.
  }
  \item{lam2PriorPar}{
    A two-element numeric vector giving the prior \eqn{\chi}{chi} and
    \eqn{\psi}{psi} parameters, respectively, for the ridge penalty
    parameter's, \eqn{\lambda_2}{lambda2}, generalized inverse gaussian
    prior.
  }
  \item{missCode}{
    An optional integer-valued code used to flag the missing data in
    \code{data}. Should take a value that cannot naturally occur in
    \code{data}. Not needed when the missing data are coded as \code{NA}.
  }
  \item{verbose}{
    A logical switch: Should verbose output be printed to stdout?\cr
    Defaults to \code{verbose = TRUE}.
  }
  \item{seed}{
    An integer-valued seed for the pseudo-random number generator. When
    \code{seed = NULL} \R's default PRNG and seed are left alone.
  }
  \item{control}{
    A list of control parameters for the Gibbs sampler and penalty parameter
    optimization (see Details for more information).
  }
}
\details{
  The \code{sampleSizes} argument must be a three-element list when
  \code{doMcem = TRUE}. In this case, the list must contain three
  two-element numeric (integer) vectors giving the number of MCMC draws
  to discard as burn-in and to retain, respectively, during the MCEM
  approximation, tuning, and sampling phases.

  When \code{doMcem = FALSE}, only a single Gibbs sample is drawn, so
  the\code{sampleSizes} argument must be a   two-element numeric
  (integer) vector giving the number of MCMC draws   to discard as
  burn-in and the number to retain in the final Gibbs sample.
  
  \code{control} is a list containing the following named elements:
  \describe{
    \item{convThresh:}{
      The R-Hat value used to judge convergence. R-Hat values <
      \code{convThresh} arising from the final, retained Gibbs sample
      will trigger a warning.\cr
      Defaults to \code{convThres = 1.1}.
    }
    \item{lambda1Starts:}{
      An optional numeric vector giving starting values for
      the LASSO penalty parameter, \eqn{\lambda_1}{lambda1}. Values are recycled
      to populate a vector with \code{size = length(targetVars)}.\cr
      Defaults to \code{rep(0.5, length(targetVars))}.
    }
    \item{lambda2Starts:}{
      An optional numeric vector giving starting values for
      the ridge penalty parameter, \eqn{\lambda_2}{lambda2}. Values are recycled
      to populate a vector with \code{size = length(targetVars)}.\cr
      Defaults to \code{rep(0.1 * nPreds, length(targetVars))}, where
      \code{nPreds} is the number of predictors in the model.
    }
    \item{usePcStarts:}{
      A logical switch: Use the starting values for \eqn{\lambda_1}{lambda1}
      suggested by Park and Casella (2008)?\cr
      Defaults to \code{usePcStarts = FALSE}.
    }
    \item{smoothingWindow:}{
      An integer giving the number of approximation phase
      \eqn{\Lambda}{Lambda} values to average over to get the starting
      \eqn{\Lambda's}{Lambdas} for the MCEM tuning phase. Setting
      \code{smoothingWindow} > 1 can facilitate convergence of the MCEM tuning phase
      when burn-in \eqn{\Lambda}{Lambda} estimates are very
      noisy. Ignored when \code{doMcem = FALSE}.\cr
      Defaults to \code{smoothingWindow = min(10, ceiling(nApprox / 10))} where
      \code{nApprox} is the number of MCEM approximation iterations. 
    }
    \item{center:}{
      A logical switch: Should the data be centered before
      estimating the imputation model? When \code{center = TRUE} the data
      centers are added back to the imputed data before the function
      returns.\cr
      Defaults to \code{center = TRUE}.
    }
    \item{scale:}{
      A logical switch: Should the predictor data be scaled to have
      unit variance before estimating the imputation model? When \code{scale =
	TRUE} imputed data are reverted to their original scaling before the
      function returns.\cr
      Defaults to \code{scale = TRUE}.
    }
    \item{adaptScales:}{
      A logical switch: Should the target variables' scales be
      actively updated as part of imputation model estimation?\cr
      Defaults to \code{adaptScales = TRUE}.
    }
    \item{simpleIntercept:}{
      A logical switch: When \code{simpleIntercept =
	TRUE}, the mean of each intercept's posterior distribution is taken as
      \eqn{\bar{y}}{\code{mean(y)}}, otherwise it equals \eqn{\bar{y} -
	\bar{\mathbf{X}}\hat{\beta}}{\code{mean(y) - colMeans(X) \%*\% Beta}}.
    }
    \item{minPredCor:}{
      The minimum correlation used by \code{mice::quickpred}
      when temporarily filling missing data before scaling or when filling
      missing data on covariates.\cr
      Defaults to \code{minPredCor = 0.3}.
    }
    \item{miceIters:}{
      The number of iterations used by the \pkg{mice} package when temporarily
      filling missing data before scaling or filling missing data on
      covariates.\cr
      Defaults to \code{miceIters = 10}.
    }
    \item{miceRidge:}{
      The ridge penalty used by the \pkg{mice} package when temporarily
      filling missing data before scaling or filling missing data on
      covariates.\cr
      Defaults to \code{miceRidge = 1e-4}.
    }
    \item{miceMethod:}{
      The elementary imputation method used by the \pkg{mice} package when
      temporarily filling missing data before scaling or filling missing data on
      covariates.\cr
      Defaults to \code{miceMethod = "pmm"}.
    }
    \item{fimlStarts:}{
      A logical switch: Should the model moments from a saturated FIML model be
      used to scale the target variables? When \code{fimlStarts = TRUE}, the
      saturated model is estimated using \pkg{lavaan}.\cr
      Defaults to \code{fimlStarts = FALSE}.
    }
    \item{optTraceLevel:}{
      A non-negative integer passed to the \pkg{optimx} \code{trace}
      argument. See \pkg{optimx} documentation for details.\cr
      Defaults to \code{optTraceLevel = 0}.
    }
    \item{optCheckKkt:}{
      A logical flag: Should the Kuhn, Jarush, Tucker optimality conditions be
      checked when optimizing the penalty parameters?\cr
      Defaults to \code{optCheckKkt = TRUE}.
    }
    \item{optMethod:}{
      A character vector giving the optimization method(s) used by \pkg{optimx}
      to estimate the penalty parameters. Possible options are
      \dQuote{Nelder-Mead}, \dQuote{BFGS}, \dQuote{CG}, \dQuote{L-BFGS-B},
      \dQuote{nlm}, \dQuote{nlminb}, \dQuote{spg}, \dQuote{ucminf},
      \dQuote{newuoa}, \dQuote{bobyqa}, \dQuote{nmkb}, \dQuote{hjkb},
      \dQuote{Rcgmin}, or \dQuote{Rvmmin}. When \code{length(optMethod) > 1},
      \pkg{optimx}'s follow-on optimization is employed. See the \pkg{optimx}
      documentation for details.\cr
      Defaults to \code{optMethod = "L-BFGS-B"}.
    }
    \item{optBoundLambda:}{
      A logical switch: Should the penalty parameters be bounded below
      by zero?\cr
      Defaults to \code{optBoundLambda = TRUE}.
    }
    \item{checkConv:}{
      A logical switch: Should the stationary Gibbs samples be checked
      for convergence?\cr
      Defaults to \code{checkConv = TRUE}.
    }
  }
}
\value{
  A reference class object with class \code{MibrrFit}. This object
  contains a great deal of metadata and various pieces of output. Key
  output can be accessed via the convenience function
  \code{\link{getParams}}, and the \code{\link{inspect}} function. See
  documentation for \code{\link{MibbrFit}} for more details.
}
\references{
  Lang, K. M. (2015) \emph{MIBEN: Multiple imputation with the
    Bayesian elastic net} (Unpublished doctoral dissertation). University of
  Kansas.
  
  Li, Q. and Lin, N. (2010) The Bayesian Elastic Net. \emph{Bayesian
    Analysis}, \bold{5}(1), 151--170.
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
  \code{\link{bl}}, \code{\link[optimx]{optimx}}, \code{\link{getParams}}, \code{\link{inspect}}
}
\examples{
data(mibrrExampleData)

## MCEM estimation:
benOut <- ben(data       = mibrrExampleData,
              y          = "y",
              X          = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
              iterations = c(30, 10)
              )

## Fully Bayesian estimation:
benOut <- ben(data         = mibrrExampleData,
              y            = "y",
              X            = setdiff(colnames(mibrrExampleData), c("y", "idNum")),
              doMcem       = FALSE,
              sampleSizes  = c(500, 500),
              lam1PriorPar = c(1.0, 0.1),
              lam2PriorPar = c(1.0, 0.1)
              )
}