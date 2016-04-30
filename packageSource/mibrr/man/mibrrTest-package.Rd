\name{mibrrTest-package}
\alias{mibrrTest-package}
\alias{mibrrTest}
\docType{package}
\title{
   Multiple Imputation with Bayesian Regularized Regression (TESTING)
}
\description{
  This package will implement a MICE-type multiple imputation method
  that uses Bayesian regularized regression models as the elementary
  imputation methods.
}
\details{
  \tabular{ll}{
    Package: \tab mibrrTest\cr
    Type: \tab Package\cr
    Version: \tab 0.1.0\cr
    Date: \tab 2015-02-23\cr
    License: \tab GPL-3\cr
  }
}
\author{
  Kyle M. Lang
  
  Maintainer: Kyle M. Lang <kylelang@ku.edu>
}
\references{
  Lang, K. M. (2015). \emph{MIBEN: Multiple imputation with the
    Bayesian elastic net} (Unpublished doctoral
  dissertation). University of Kansas.
  
  Li, Q., & Lin, N. (2010). The Bayesian Elastic Net. \emph{Bayesian
    Analysis, 5}(1), 151--170.
  
  Park, T., & Casella, G. (2008). The Bayesian Lasso. \emph{Journal of
    the American Statistical Association, 103}, 681--686.
  
  Zhao, Y., & Long, Q. (2013). Multiple imputation in the presence of
  high-dimensional data. \emph{Statistical Methods in Medical Research,
    0}(0), 1--15.
}
\keyword{ package }
\examples{
data(mibrrExampleData)

miben(rawData = mibrrExampleData,
      targetVars = c("y", paste0("x", c(1 : 3))),
      ignoreVars = "idNum")

mibl(rawData = mibrrExampleData,
     targetVars = c("y", paste0("x", c(1 : 3))),
     ignoreVars = "idNum")
}
