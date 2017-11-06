\name{mibl}
\alias{mibl}
\title{
Multiple Imputation with the Bayesian LASSO
}
\description{
  This function will implement a MICE-based multiple imputation scheme that
  employs the Bayesian LASSO as its elementary imputation method.
}
\usage{
mibl(data,
     nImps,
     targetVars     = NULL,
     ignoreVars     = NULL,
     iterations     = c(100, 10),
     sampleSizes    = list(rep(25, 2), rep(250, 2), rep(500, 2)),
     missCode       = NULL,
     returnConvInfo = TRUE,
     returnParams   = FALSE,
     verbose        = TRUE,
     seed           = NULL,
     control        = list()
     )
}
\arguments{
  \item{data}{
    An incomplete, numeric data matrix or data frame for which to create the
    imputations.
  }
  \item{nImps}{
    An integer giving the number of imputations to create.
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
    during the MCEM approximation and tuning phases, respectively.\cr
    Defaults to \code{iterations = c(100, 10)}.
  }
  \item{sampleSizes}{
    A list containing three two-element numeric vectors giving the number of
    MCMC draws discarded as burn-in and retained, respectively, during the MCEM
    approximation, tuning, and sampling phases.\cr
    Defaults to \code{sampleSizes = list(rep(25, 2), rep(250, 2), rep(500, 2))}. 
  }
  \item{missCode}{
    An optional integer-valued code used to flag the missing data in
    \code{data}. Should take a value that cannot naturally occur in
    \code{data}. Not needed when the missing data are coded as \code{NA}.
  }
  \item{returnConvInfo}{
    A logical switch: Should convergence information for the imputation model
    (i.e., history of the optimized penalty parameters and R-Hat values for
    final parameter estimates) be returned?\cr
    Defaults to \code{returnConvInfo = TRUE}. 
  }
  \item{returnParams}{
    A logical switch: Should the final Gibbs samples of the imputation model's
    parameters be returned?\cr
    Defaults to \code{returnParams = FALSE}.
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
  \code{control} is a list containing the following named elements:
  \describe{
    \item{convThresh:}{
      The R-Hat value used to judge convergence. R-Hat
      values < \code{convThresh} arising during the MCEM tuning phase will
      trigger a warning.\cr
      Defaults to \code{convThres = 1.1}.
    }
    \item{lambda1Starts:}{
      An optional numeric vector giving starting values for
      the LASSO penalty parameter, \eqn{\lambda_1}{lambda1}. Values are recycled
      to populate a vector with \code{size = length(targetVars)}.\cr
      Defaults to \code{rep(0.5, length(targetVars))}.
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
      when burn-in \eqn{\Lambda}{Lambda} estimates are very noisy.\cr
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
    \item{preserveStructure:}{
      A logical switch: Should the data columns be returned in the same order
      as submitted?\cr
      Defaults to \code{preserveStructure = TRUE}.
    }
  }
}
\value{
  A list containing \code{nImps} imputed versions of \code{data} and
  whatever additional output is requested via \code{returnConvInfo} and
  \code{returnParams}. 
}

\references{
  Lang, K. M. (2015) \emph{MIBEN: Multiple imputation with the
    Bayesian elastic net} (Unpublished doctoral dissertation). University of
  Kansas.
  
  Park, T. and Casella, G. (2008) The Bayesian Lasso. \emph{Journal of
    the American Statistical Association}, \bold{103}, 681--686.
  
  Zhao, Y., and Long, Q. (2013) Multiple imputation in the presence of
  high-dimensional data. \emph{Statistical Methods in Medical Research},
  \bold{0}(0), 1--15.
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
\code{\link{miben}}
}
\examples{
data(mibrrExampleData)

miblOut <- mibl(data       = mibrrExampleData,
                nImps      = 100,
                iterations = c(50, 10),
                targetVars = c("y", paste0("x", c(1 : 3))),
                ignoreVars = "idNum")
}