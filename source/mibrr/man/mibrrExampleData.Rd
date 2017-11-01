\name{mibrrExampleData}
\alias{mibrrExampleData}
\docType{data}
\title{
Example Dataset for the MIBRR package.
}
\description{
Toy data generated as in Experiment 1 of Lang (2015).
}
\usage{data("mibrrExampleData")}
\format{
  A data frame with 200 observations on the following 17 variables.
  \describe{
    \item{idNum:}{
      \dQuote{Participant} ID Number.
    }
    \item{y:}{
      Outcome Variable. Contains 20\% MAR missingness.
    }
    \item{x1--x3:}{
      Substantive predictors. Contain 20\% MAR missingness.
    }
    \item{z1--z12:}{
      Exogenous auxiliary variables. Contain 10\% MCAR missingness. 
    }
  }
}
\details{
These data are only simulated toy data; they have no true meaning. This
is one of the datasets generated as part of the Monte Carlo simulation
study used to conduct Experiment 1 of my dissertation. The missingness
on \{y, X\} is caused by a linear combination of two randomly selected
elements of \{Z\}, and half of the elements in \{Z\} have no association
with \{y, X\} (see Lang, 2015, for more details). 
}
\source{
  Lang, K. M. (2015) \emph{MIBEN: Multiple imputation with the Bayesian
    elastic net} (Unpublished doctoral dissertation). University of
  Kansas.
}
\keyword{datasets}
