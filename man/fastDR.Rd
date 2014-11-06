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
       verbose=FALSE)
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

\item{verbose}{If TRUE, the program will print out progress and performance
indicators}

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

\code{fastDR} uses a generalized boosted model to estimate a propensity score
model (McCaffrey, Ridgeway, and Morral, 2004). It uses those propensity scores
as weights in a weighted generalized linear model to produce a doubly robust
estimate of the treatment effect.

To reuse the propensity score weights in other analyses, make sure that the code aligns the right weight with the right case. Use the key in the data and the key attached to the weights. For example, \code{mydata$psw <- myfastDR$w[match(mydata$observationID,names(myfastDR$w))]}

The effective sample size is the number of independent cases that would yield the same precision as the given weighted sample. It is computed as
\deqn{\frac{(\sum_i (1-t_i)w_i)^2}{\sum_i (1-t_i)w_i^2}}

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
\item{effects}{a table showing the results. See the details section}
\item{y.dist}{contains the value of \code{y.dist} used in the call to \code{fastDR}}
\item{shrinkage}{the shrinkage parameter used in the propensity score model}


}

\references{

D. McCaffrey, G. Ridgeway, and A. Morral (2004). \dQuote{Propensity score
estimation with boosted regression for evaluating adolescent substance abuse
treatment,} \emph{Psychological Methods} 9(4):403-425.

J.H. Friedman (2001). \dQuote{Greedy Function Approximation: A Gradient
Boosting Machine,} \emph{Annals of Statistics} 29(5):1189-1232.

\url{http://sites.google.com/site/gregridgeway}

}


\author{Greg Ridgeway \email{gridge@sas.upenn.edu}}

\seealso{ \code{\link{gbm}} }

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
              verbose=TRUE)             # print out detailed progress

# show balance table
round(dr1$balance.tab,3)

# compute odds-ratio
print(dr1,type="outcome", model="dr")

}

\keyword{models}

