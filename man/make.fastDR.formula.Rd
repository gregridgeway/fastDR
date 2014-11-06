\name{make.fastDR.formula}
\alias{make.fastDR.formula}

\title{
Make formulas appropriate for using in fastDR
}
\description{
\code{make.fastDR.formula} is a function that helps the user write appropriate formulas for use in fastDR. Since
}
\usage{
make.fastDR.formula(y.vars,
                    t.var,
                    x.vars=".",
                    weights.var=NULL,
                    key.var,
                    data=NULL)
}

\arguments{
   \item{y.vars}{a vector of strings indicating the names of the outcome variables (e.g. \code{c("y1","y2")})}
   \item{t.var}{a string giving the name of the treatment variable (e.g. \code{"treat"})}
   \item{x.vars}{a vector of strings giving the names of the covariates (e.g. \code{c("X1","X2")}). \code{x.vars} may also take the value \code{"."}, which will select all variables in \code{data} not otherwise given in the other components}
   \item{weights.var}{a string giving the name of the observation weights (e.g. \code{"samp.weights"})}
   \item{key.var}{a string giving the name of the variable containing the observation IDs (e.g. \code{"caseID"})}
   \item{data}{a data frame that will be sent to \code{fastDR}}
}

\value{
Returns a list with appropriately formatted R formulas ready to be submitted to \code{fastDR}
}

\author{Greg Ridgeway \email{gridge@sas.upenn.edu}}

\seealso{ \code{\link{fastDR}} }

\examples{
# NHANES example from survey package
data(nhanes)

# add a unique ID to each row
nhanes$observationID <- 1:nrow(nhanes)
# recode the "treatment" (male) to a 0/1 indicator
nhanes$male <- as.numeric(nhanes$RIAGENDR==1)
# create a second random outcome
nhanes$Y2 <- rnorm(nrow(nhanes),0,1)
# drop unused variables
j <- which((names(nhanes) == "SDMVPSU") |
           (names(nhanes) == "SDMVSTRA"))
nhanes <- nhanes[,-j]

my.forms <- make.fastDR.formula(y.vars=c("HI_CHOL","Y2"),
                                t.var="male",
                                x.vars=".",
                                weights.var="WTMEC2YR",
                                key.var="observationID",
                                data=nhanes)

}

\keyword{ models }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
