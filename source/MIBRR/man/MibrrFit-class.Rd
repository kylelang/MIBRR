\name{MibrrFit-class}
\Rdversion{1.1}
\docType{class}
\alias{MibrrFit-class}
\alias{MibrrFit}

\title{Class \code{"MibrrFit"}}
\description{
  The \code{MibrrFit} class is a reference class object that acts as the
  general return object for the \code{MIBRR} package.  
}
\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.
}
\author{
  Kyle M. Lang
}
\note{
  \code{MibrrFit} is a mutable \emph{reference class} object, so
  functions in the \code{MIBRR} package that take a \code{MibrrFit}
  object as input will, generally, alter the state of that input object.
  
  If you need to preserve a snapshot of a \code{MibrrFit} object, use
  the \code{copy} member function (e.g., \code{snapshot <-
    mibrrFit$copy()}).
}
\seealso{
\code{\link{getImpData}}, \code{\link{getParams}}, \code{\link{getField}},
\code{\link{postPredict}}
}
\section{Warning}{
  This function is in a highly unstable \emph{alpha} level of development. Please
  anticipate frequent---and dramatic---changes to the functionality and user
  interface.
  
  You have been granted access to this package for evaluation purposes, only. This
  function is \strong{absolutely not} ready for use in real-world analyses!
}
\examples{
showClass("MibrrFit")
}
\keyword{classes}
\section{Methods}{
  The \code{MibrrFit} class contains many private methods that are not
  meant to be directly accessed by end users. Adventurous users can
  refer to the source code of the class' definition for details.

  Any field in the \code{MibrrFit} object can be accessed via the
  \code{\link{getField}} function. The most interesting fields are
  described in the \code{\link{getField}} function's documentation.

  After running \code{miben} or \code{mibl}, a set of imputed datasets
  can be generated with the \code{\link{getImpData}} function.

  After running any of the primary modeling functions (i.e.,
  \code{miben}, \code{mibl}, \code{ben}, \code{bl}), the model
  parameters' posterior Gibbs samples can be extracted via the
  \code{\link{getParams}} function and predictions can be generated via
  the \code{\link{postPredict}} function.
}
\section{Fields}{
  The \code{MibrrFit} class contains the following fields (which are not
  intended to be directly access by users). Interesting data, results,
  and metadata can be extracted from a fitted \code{MibrrFit} object via
  the functions described in the preceding \dQuote{methods} section.
  
  \describe{
    \item{\code{data}:}{Object of class \code{data.frame} ~~ }
    \item{\code{targetVars}:}{Object of class \code{character} ~~ }
    \item{\code{ignoreVars}:}{Object of class \code{character} ~~ }
    \item{\code{iterations}:}{Object of class \code{integer} ~~ }
    \item{\code{sampleSizes}:}{Object of class \code{list} ~~ }
    \item{\code{missCode}:}{Object of class \code{integer} ~~ }
    \item{\code{seed}:}{Object of class \code{list} ~~ }
    \item{\code{doImp}:}{Object of class \code{logical} ~~ }
    \item{\code{doMcem}:}{Object of class \code{logical} ~~ }
    \item{\code{doBl}:}{Object of class \code{logical} ~~ }
    \item{\code{checkConv}:}{Object of class \code{logical} ~~ }
    \item{\code{verbose}:}{Object of class \code{logical} ~~ }
    \item{\code{convThresh}:}{Object of class \code{numeric} ~~ }
    \item{\code{lambda1Starts}:}{Object of class \code{numeric} ~~ }
    \item{\code{lambda2Starts}:}{Object of class \code{numeric} ~~ }
    \item{\code{l1Pars}:}{Object of class \code{numeric} ~~ }
    \item{\code{l2Pars}:}{Object of class \code{numeric} ~~ }
    \item{\code{usePcStarts}:}{Object of class \code{logical} ~~ }
    \item{\code{smoothingWindow}:}{Object of class \code{integer} ~~ }
    \item{\code{center}:}{Object of class \code{logical} ~~ }
    \item{\code{scale}:}{Object of class \code{logical} ~~ }
    \item{\code{adaptScales}:}{Object of class \code{logical} ~~ }
    \item{\code{minPredCor}:}{Object of class \code{numeric} ~~ }
    \item{\code{miceIters}:}{Object of class \code{integer} ~~ }
    \item{\code{miceRidge}:}{Object of class \code{numeric} ~~ }
    \item{\code{miceMethod}:}{Object of class \code{character} ~~ }
    \item{\code{fimlStarts}:}{Object of class \code{logical} ~~ }
    \item{\code{preserveStructure}:}{Object of class \code{logical} ~~ }
    \item{\code{optTraceLevel}:}{Object of class \code{integer} ~~ }
    \item{\code{optCheckKkt}:}{Object of class \code{logical} ~~ }
    \item{\code{optMethod}:}{Object of class \code{character} ~~ }
    \item{\code{optBoundLambda}:}{Object of class \code{logical} ~~ }
    \item{\code{dataMeans}:}{Object of class \code{numeric} ~~ }
    \item{\code{dataScales}:}{Object of class \code{numeric} ~~ }
    \item{\code{gibbsOut}:}{Object of class \code{list} ~~ }
    \item{\code{ignoredColumns}:}{Object of class \code{data.frame} ~~ }
    \item{\code{rawNames}:}{Object of class \code{character} ~~ }
    \item{\code{impRowsPool}:}{Object of class \code{integer} ~~ }
    \item{\code{missList}:}{Object of class \code{list} ~~ }
    \item{\code{nChains}:}{Object of class \code{integer} ~~ }
    \item{\code{rHats}:}{Object of class \code{list} ~~ }
    \item{\code{lambdaMat}:}{Object of class \code{matrix} ~~ }
    \item{\code{lambdaHistory}:}{Object of class \code{list} ~~ }
    \item{\code{betaStarts}:}{Object of class \code{matrix} ~~ }
    \item{\code{tauStarts}:}{Object of class \code{matrix} ~~ }
    \item{\code{sigmaStarts}:}{Object of class \code{numeric} ~~ }
    \item{\code{userMissCode}:}{Object of class \code{logical} ~~ }
    \item{\code{missCounts}:}{Object of class \code{integer} ~~ }
    \item{\code{nTargets}:}{Object of class \code{integer} ~~ }
    \item{\code{nVar}:}{Object of class \code{integer} ~~ }
    \item{\code{nPreds}:}{Object of class \code{integer} ~~ }
    \item{\code{nObs}:}{Object of class \code{integer} ~~ }
    \item{\code{totalIters}:}{Object of class \code{integer} ~~ }
  }
}

