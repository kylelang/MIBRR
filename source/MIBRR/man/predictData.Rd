\name{predictData}
\alias{predictData}
\docType{data}
\title{
Example Dataset to Demonstrate Prediction with the MIBRR Package.
}
\description{
These data can be used to demonstrate the prediction capabilities of
\code{ben}, \code{bl}, \code{miben}, and \code{mibl}.
}
\usage{data("predictData")}
\format{
  A three-element list containing the data frames described below.
  
  All three data frames contain the same 11 variables. Scores
  on five personality dimensions (i.e., agreeableness,
  conscientiousness, extroversion, neuroticism, and openness to
  experience). Participants' ages, genders, and a set of four dummy
  codes indicating their level of educational achievement.
  \describe{
    \item{incomplete:}{
      An incomplete dataset with 500 observations of 11 variables. The
      five personality scores in this dataset are each subject of
      approximately 20\% missing data.
    }
    \item{test:}{
      A testing dataset with 50 observations of 11 variables.
    }
    \item{train:}{
      A training dataset with 500 observations of 11 variables.
    }
  }
}
\details{
  These data are derived from the \code{bfi} data provided by the
  \pkg{psych} package. I scored the five personality constructs, dummy
  coded the covariates, and sub-sampled 500 observations. The incomplete
  dataset had 20\% MAR missing data imposed on the personality
  constructs by defining the reponse propensity as a function of the
  fully observed covariates.
}
\source{
  Revelle, W. (2017) psych: Procedures for Personality and
  Psychological Research, Northwestern University, Evanston, Illinois,
  USA, https://CRAN.R-project.org/package=psych Version = 1.7.8.
}
\keyword{datasets}

