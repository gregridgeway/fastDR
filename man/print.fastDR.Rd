\name{print.fastDR}
\alias{print.fastDR}

\title{ Print model summary }
\description{
  Display basic information about a \code{fastDR} object.
}
\usage{
\method{print}{fastDR}(x,type="outcome",model="dr",...)
}
\arguments{
  \item{x}{ an object of class \code{fastDR}}
  \item{type}{a string taking value \dQuote{outcome} or \dQuote{complete}}
  \item{model}{a string taking value \dQuote{un}, \dQuote{ps}, or \dQuote{dr}}
  \item{\dots}{ arguments passed to \code{print.default}}
}
\details{
If \code{type="outcome"} then \code{print.fastDR()} will output the treatment effect estimate. Which treatment effect estimate depends on the value of \code{model}. Presumably, \code{"dr"} will be needed, but an unadjusted and propensity score weighted option are available. The results will be transformed to a scale appropriate for the \code{y.dist} specified when the user called \code{fastDR}. \code{print.fastDR()} will produce odds ratios for binomial outcomes, rate ratios for Poisson outcomes, and differences for all others.

If \code{type="complete"} then \code{print.fastDR()} will print the entire effects results.

See \code{\link{fastDR}} for an example of its use.
}


\author{ Greg Ridgeway \email{gridge@sas.upenn.edu}}

\seealso{ \code{\link{fastDR}} }

\keyword{models}
