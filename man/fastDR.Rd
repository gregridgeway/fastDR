\name{fastDR}
\alias{fastDR}
\title{Fast Doubly Robust Estimation}
\description{Doubly robust treatment effect estimation with non-parametric propensity score estimates}
\usage{
fastDR(form.list,
       data,
       y.dist="gaussian",
       n.trees=3000,
       interaction.depth=3,
       shrinkage=0.01,
       verbose=FALSE,
       ps.only=FALSE,
       keepGLM=TRUE,
       smooth.lm=0,
       par_details=gbmParallel(1,1024))
}

\arguments{

\item{form.list}{a list of formulas giving the components of the model, outcomes, treatment indicator, covariates, observation weights, and observation keys. Consider using \code{\link{make.fastDR.formula}} to assist in creating \code{form.list}}

\item{data}{a data frame containing the variables in the model}

\item{y.dist}{the distribution assumed for the outcome. These should be
specified as character strings (e.g. "gaussian", "quasibinomial",
"quasipoisson"). Using "binomial" will cause inconsequential warning messages
since the model uses non-integer propensity score weights. \dQuote{quasibinomial}
avoids this. If \code{y.form} has multiple outcomes than \code{y.dist} should
be a vector of strings the same length as the number of outcomes in
\code{y.dist}}

\item{n.trees}{the total number of trees to fit for the generalized boosted
model. This is equivalent to the number of iterations and the number of basis
functions in the additive expansion. The default is most likely appropriate}

\item{interaction.depth}{The maximum depth of variable interactions. 1 implies
an additive model, 2 implies a model with up to 2-way interactions, etc.}

\item{shrinkage}{a shrinkage parameter applied to each tree in the generalized
boosted model expansion. Also known as the learning rate or step-size
reduction. The default is most likely appropriate}

\item{verbose}{If TRUE, the fastDR will print out progress and performance
indicators}

\item{ps.only}{If TRUE, the fastDR will skip the DR step and just return the propensity scores. This is useful if you need to do non-standard changes to the regression model, like eliminating some covariates}

\item{keepGLM}{if TRUE, keep a copy of the \code{glm} objects produced in the outcome model step. These models can become large with multiple outcomes and large datasets}

\item{smooth.lm}{a numeric value for penalizing the size of the covariates in the DR step}

\item{par_details}{a list containing the number of threads and array_chunk_size to be passed to \code{gbmt} for running several gbm computations in parallel (e.g. gradient calculation, split selection). Most easily set using \code{gbmParallel}}

}


\details{
The form.list has the following components.

\code{y.form} is a formula specifying the outcomes (e.g. \code{~y}). Multiple
outcomes may be listed (e.g. \code{~y1+y2}); fastDR will conduct a separate
analysis for each of them

\code{t.form} is a formula specifying the treatment indicator (e.g. ~treat). The
formula can include only one binary treatment coded as 0/1

\code{x.form} is a formula specifying the potential confounding variables that
will be included in the propensity score model and in the outcome model (e.g.
\code{~X1+X2+X2})

\code{weights.form} is an optional formula specifying observation weights, such as sampling weights

\code{key.form} is a formula giving observation IDs (e.g. \code{~caseID}). This is required in order to make sure that after returning results, the propensity
score weights can be correctly aligned with cases by the user

\code{fastDR} only estimates the average treatment effect on the treated.

\code{fastDR} uses a generalized boosted model to estimate a propensity score model (McCaffrey, Ridgeway, and Morral, 2004). It uses those propensity scores as weights in a weighted generalized linear model to produce a doubly robust estimate of the treatment effect.

To reuse the propensity score weights in other analyses, make sure that the code aligns the right weight with the right case. Use the key in the data and the key attached to the weights. For example, \code{mydata$psw <- myfastDR$w[match(mydata$observationID,names(myfastDR$w))]}

The effective sample size is the number of independent cases that would yield the same precision as the given weighted sample. It is computed as
\deqn{\frac{(\sum_i (1-t_i)w_i)^2}{\sum_i (1-t_i)w_i^2}}

When \code{smooth.lm} is greater than 0, \code{fastDR} uses a data augmentation method to impose a penalty on the regression coefficients in the DR step. At a minimum, modest penalties provide numerical stability in the event of correlation among the covariates or certain features being highly predictive of the outcome. \code{fastDR} does not penalize the intercept or the treatment indicator. \code{fastDR} uses the data augmentation approach described in Greenland and Mansournia (2015). This approach uses ridge regression penalties (N(0,m)) for Gaussian outcomes, log-F(m,m) prior for logistic regression, and gamma(m,m) prior for Poisson regression. The appropriate size of the penalty is subjective, but values between 0.1 and 4.0 are modest penalties for general application.

If the number of treatment cases exceeds 1500, then \code{fastDR} will take a random sample of 1500 cases for computing the standard errors of the estimates of E(Y1|t=1) and E(Y1-Y0|t=1). Computing the standard error requires computing the mean of the n1*n1 covariance matrix of predicted values, which can require a lot of memory when n1 is large, but can be efficiently estimated with a sample. Note that there will be some amount of sampling variability.

}

\value{ \code{fastDR} returns a \code{fastDR.object}, a list containing

\item{ks}{the largest KS statistic in the balance table (see \code{balance.tab} below). A measure of the quality of the propensity score fit}
\item{ks.un}{the largest KS statistic in the balance table before propensity score weighting (see \code{balance.tab.un} below)}
\item{balance.tab}{the balance table showing the treatment means and propensity score weighted control means for each of the terms specified in \code{x.form}. For categorical features the table lists each level separately}
\item{balance.tab.un}{the balance table showing the treatment means and unweighted control means for each of the terms specified in \code{x.form}}
\item{w}{the propensity score weights used in the analysis. Use \code{names()} to extract the key associated with each weight to align with the original data. All cases in the original data might not necessarily have a propensity score weight (e.g. missing value for sampling weight)}
\item{p}{the estimated propensity scores. \code{p} also maintains the \code{key} in its names}
\item{n1,ESS}{the number of cases in the treatment group and the effective sample size of the control group}
\item{glm.un,glm.ps,glm.dr}{the generalized linear model fits using only sampling weights (\code{glm.un}), using propensity score weights (\code{glm.ps}), and the model producing doubly robust estimates (\code{glm.dr})}
\item{effects}{a table showing the results. \code{effects} includes results for an unadjusted (un), propensity score adjusted (ps), and doubly robust estimate (dr). For each, \code{effects} includes the estimates and standard errors for E(Y1|t=1), E(Y0|t=1), their difference, and p-value. These results are all on the scale of the response (not odds ratios or rate ratios)}
\item{y.dist}{contains the value of \code{y.dist} used in the call to \code{fastDR}}
\item{shrinkage}{the shrinkage parameter used in the propensity score model}


}

\references{

D. McCaffrey, G. Ridgeway, and A. Morral (2004). \dQuote{Propensity score
estimation with boosted regression for evaluating adolescent substance abuse
treatment,} \emph{Psychological Methods} 9(4):403-425.

J.H. Friedman (2001). \dQuote{Greedy Function Approximation: A Gradient
Boosting Machine,} \emph{Annals of Statistics} 29(5):1189-1232.

S. Greenland and M.A. Mansournia (2015). \dQuote{Penalization, bias reduction, and default priors in logistic and related categorical and survival regressions,} \emph{Statistics in Medicine} 34:3133-3143.

S. Greenland, M.A. Mansournia, and D.G. Altman (2014). \dQuote{Sparse data bias: a problem hiding in plain sight,} \emph{The BMJ} 353:i1981.

\url{http://sites.google.com/site/gregridgeway}

}


\author{Greg Ridgeway \email{gridge@upenn.edu}}

\seealso{ \code{\link{gbmt}}, \code{\link{gbmParallel}} }

\examples{
# NHANES example from survey package

data(nhanes)
# add a unique ID to each row
nhanes$observationID <- 1:nrow(nhanes)
# recode the "treatment" (male) to a 0/1 indicator
nhanes$male <- as.numeric(nhanes$RIAGENDR==1)
# make "race" a factor with descriptive labels
nhanes$race <- factor(nhanes$race,
                      levels=c(1,2,3,4),
                      labels=c("Hispanic","non-Hispanic white",
                               "non-Hispanic black","other"))

dr1 <- fastDR(list(y.form=~HI_CHOL,          # high cholesterol (over 240mg/dl)
                   t.form=~male,             # compare males to similar females
                   x.form=~race+agecat,      # potential confounders
                   weights.form=~WTMEC2YR,   # sampling weights
                   key.form=~observationID),
              data=nhanes,
              y.dist="quasibinomial",   # outcome is 0/1
              n.trees=1000,
              interaction.depth=3,
              shrinkage=0.005,
              verbose=TRUE,             # print out detailed progress
              smooth.lm=0.1,            # stabilize regression coefs
              par_details=gbmParallel(1,1024)) # just use one core

# show balance table
round(dr1$balance.tab,3)

# compute percentage difference
print(dr1, type="outcome", model="dr")

}

\keyword{models}
