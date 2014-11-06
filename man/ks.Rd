\name{ks}
\alias{ks}
\title{Weighted Kolmogorov-Smirnov Statistic}
\description{Computes a weighted Kolmogorov-Smirnov statistic to measure the difference in the marginal feature distributions between the treatment and control cases}
\usage{
ks(x, z, w)
}

\arguments{
  \item{x}{a vector of numeric measurements}
  \item{z}{a vector of 0/1 indicators indicating group membership}
  \item{w}{a vector of weights}
}

\details{
This function is used in the propensity score model building step to assess the balance between the treatment and control cases
}

\value{
Returns the Kolmogorov-Smirnov statistic, the largest difference between treatment and control groups' empirical cumulative distribution functions
}

\references{

A. Kolmogorov (1933). \dQuote{Sulla determinazione empirica di una legge di distribuzione,}. \emph{Giornale dell'Istituto Italiano degli Attuari} 4:83-91.

N. Smirnov (1948). \dQuote{Table for estimating the goodness of fit of empirical distributions,}. \emph{Annals of Mathematical Statistics} 19:279-281.

}


\author{Greg Ridgeway \email{gridge@sas.upenn.edu}}

\seealso{ \code{\link{fastDR}} }

\examples{
y     <- c(rnorm(100,0,1),rnorm(100,0.5,1))
treat <- rep(0:1,each=100)
w     <- 1/c(pnorm(y[1:100],0,1),pnorm(y[101:200],0.5,1))
ks(x=y,z=treat,w=w)
}

\keyword{univar}
