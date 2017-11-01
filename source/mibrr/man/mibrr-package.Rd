\name{mibrr-package}
\alias{mibrr-package}
\docType{package}
\title{
\packageTitle{mibrr}
}
\description{
\packageDescription{mibrr}
}
\details{
\packageIndices{mibrr}
}
\author{
\packageAuthor{mibrr}

Maintainer: \packageMaintainer{mibrr}
}
\references{
  Lang, K. M. (2015) \emph{MIBEN: Multiple imputation with the Bayesian
    elastic net} (Unpublished doctoral dissertation). University of
  Kansas.
  
  Li, Q. and Lin, N. (2010) The Bayesian Elastic Net. \emph{Bayesian
  Analysis}, \bold{5}(1), 151--170.
  
  Park, T. and Casella, G. (2008) The Bayesian Lasso. \emph{Journal of
    the American Statistical Association}, \bold{103}, 681--686.
  
  Zhao, Y., and Long, Q. (2013) Multiple imputation in the presence of
  high-dimensional data. \emph{Statistical Methods in Medical Research},
    \bold{0}(0), 1--15.
}
\examples{
data(mibrrExampleData)

miben(data       = mibrrExampleData,
      nImps      = 100,
      targetVars = c("y", paste0("x", c(1 : 3))),
      ignoreVars = "idNum")

mibl(data       = mibrrExampleData,
     nImps      = 100,
     targetVars = c("y", paste0("x", c(1 : 3))),
     ignoreVars = "idNum")
}
\keyword{ package }
