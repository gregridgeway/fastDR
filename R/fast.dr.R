# Notes:
#   - fix treatment cases weighted with sampling weights

#' @importFrom gbm3 gbm_dist gbmParallel gbmt training_params
#' @importFrom methods is
#' @importFrom stats coef cor delete.response formula median model.frame
#'   model.matrix na.pass predict reformulate terms update vcov
#'   weighted.mean
#' @importFrom survey svydesign svyglm svymean
"_PACKAGE"

.onAttach <- function(libname, pkgname)
{
   packageStartupMessage(paste("Loaded fastDR",
                               utils::packageDescription("fastDR")$Version))
   
   # double check survey package version... critical
   if(utils::packageDescription("survey")$Version < "4.1")
   {
      packageStartupMessage('WARNING: fastDR requires version 4.1 or later of the survey package. Earlier versions of the survey package did not have a rescale=FALSE option or had a bug when rescale=FALSE that did not copy the weights into the dataset. install.packages("survey") to obtain the latest version.')
   }
}

#' Make formulas appropriate for using in fastDR
#'
#' `make.fastDR.formula()` helps the user write appropriate formulas for use in
#' [fastDR()].
#'
#' @param y.vars A vector of strings indicating the names of the outcome
#'   variables, for example `c("y1", "y2")`.
#' @param t.var A string giving the name of the treatment variable, for example
#'   `"treat"`.
#' @param x.vars A vector of strings giving the names of the covariates, for
#'   example `c("X1", "X2")`. `x.vars` may also take the value `"."`, which
#'   will select all variables in `data` not otherwise given in the other
#'   components.
#' @param weights.var A string giving the name of the observation weights, for
#'   example `"samp.weights"`.
#' @param key.var A string giving the name of the variable containing the
#'   observation IDs, for example `"caseID"`.
#' @param data A data frame that will be sent to [fastDR()].
#'
#' @return A list with appropriately formatted R formulas ready to be submitted
#'   to [fastDR()].
#'
#' @author Greg Ridgeway \email{gridge@upenn.edu}
#' @seealso [fastDR()]
#'
#' @examples
#' # NHANES example from survey package
#' data(nhanes)
#'
#' # add a unique ID to each row
#' nhanes$observationID <- 1:nrow(nhanes)
#' # recode the "treatment" (male) to a 0/1 indicator
#' nhanes$male <- as.numeric(nhanes$RIAGENDR == 1)
#' # create a second random outcome
#' nhanes$Y2 <- rnorm(nrow(nhanes), 0, 1)
#' # drop unused variables
#' j <- which((names(nhanes) == "SDMVPSU") |
#'            (names(nhanes) == "SDMVSTRA"))
#' nhanes <- nhanes[, -j]
#'
#' my.forms <- make.fastDR.formula(y.vars = c("HI_CHOL", "Y2"),
#'                                 t.var = "male",
#'                                 x.vars = ".",
#'                                 weights.var = "WTMEC2YR",
#'                                 key.var = "observationID",
#'                                 data = nhanes)
#'
#' @keywords models
#' @export
make.fastDR.formula <-
   function(y.vars,
            t.var,
            x.vars=".",
            weights.var=NULL,
            key.var,
            data=NULL)
{
   if(!is.null(data))
   {
      if(!all(y.vars %in% names(data)))
         warning("Not all y.vars appear in data")
      if(!(t.var %in% names(data)))
         warning("t.var does not appear in data")
      if(!is.null(weights.var) && !(weights.var %in% names(data)))
         warning("weights.var does not appear in data")
      if(!(key.var %in% names(data)))
         warning("key.var does not appear in data")
      if(identical(x.vars, "."))
         x.vars <- setdiff(names(data),c(y.vars,t.var,weights.var,key.var))
   }
   if(is.null(data) && identical(x.vars, "."))
      stop("If x.vars given as '.' then data must be given")

   return(list(
      y.form      =formula(paste("~",paste(y.vars,collapse="+"),sep="")),
      t.form      =formula(paste("~",t.var,sep="")),
      x.form      =formula(paste("~",paste(x.vars,collapse="+"),sep="")),
      weights.form=if(is.null(weights.var)) NULL else formula(paste("~",weights.var,sep="")),
      key.form    =formula(paste("~",key.var,sep=""))))
}

#' Print a fastDR object
#'
#' Display basic information about a `fastDR` object.
#'
#' @param x An object of class `fastDR`.
#' @param type A string taking value `"outcome"` or `"complete"`.
#' @param model A string taking value `"un"`, `"ps"`, or `"dr"`.
#' @param ... Arguments passed to `print.default`.
#'
#' @details
#' If `type = "outcome"` then `print.fastDR()` will output the treatment effect
#' estimate. Which treatment effect estimate depends on the value of `model`.
#' Presumably, `"dr"` will be needed, but an unadjusted and propensity score
#' weighted option are available. Results are reported on the response scale,
#' not as odds ratios or rate ratios. `print.fastDR()` will produce percentage
#' point differences for binomial outcomes, rate differences for Poisson
#' outcomes, and differences for all others.
#'
#' If `type = "complete"` then `print.fastDR()` will print the entire effects
#' results.
#'
#' See [fastDR()] for an example of its use.
#'
#' @return `print.fastDR()` is called for its side effect of printing a summary
#'   of treatment effect estimates. It returns `NULL` invisibly.
#'
#' @author Greg Ridgeway \email{gridge@upenn.edu}
#' @seealso [fastDR()]
#'
#' @keywords models
#' @method print fastDR
#' @export
print.fastDR <- function(x, type="outcome", model="dr",... )
{
   if(!inherits(x, "fastDR"))
   {
      stop("object must be a fastDR object, typically produced from a call to fastDR")
   }

   if(!(type %in% c("outcome","complete")))
      stop("type parameter must be one of 'outcome' (default) or 'complete'")

   if(is.null(x$effects))
   {
      stop("The fitted fastDR model has no estimated effects. Most likely the call to fastDR had ps.only=TRUE, which only estimates propensity scores and does not estimate treatment effects")
   } else if(type=="outcome")
   {
      if(!(model %in% c("un","ps","dr")))
         stop("model parameter must be one of 'un', 'ps', or 'dr' (default)")
      if(model=="un")
         cat("Unadjusted results\n")
      else if(model=="ps")
         cat("Propensity score weighting results\n")
      else if(model=="dr")
         cat("Doubly robust results\n")
      cat("Estimand: ", x$estimand, "\n")
      for(i in 1:length(x$effects))
      {
         cat(names(x$effects)[i],":\n",sep="")
         temp <- x$effects[[i]]
         temp <- temp[model,"TE"] + c(0,-2,2)*temp[model,"se.TE"]
         temp <- signif(temp,3)
         if(x$y.dist[i] %in% c("binomial","quasibinomial"))
         {
            cat("percentage point difference (95% CI): ",100*temp[1],
                " (",100*temp[2],",",100*temp[3],")\n",sep="")
         }
         else if(x$y.dist[i] %in% c("poisson","quasipoisson"))
         {
            cat("rate difference (95% CI): ",temp[1],
                " (",temp[2],",",temp[3],")\n",sep="")
         }
         else
         {
            cat("effect (95% CI): ",temp[1],
                " (",temp[2],",",temp[3],")\n",sep="")
         }
      }
   }

   if(type=="complete")
   {
      print(x$effects)
   }
}

#' Summarize a fastDR object
#'
#' Display treatment effect estimates and balance diagnostics for a `fastDR`
#' object.
#'
#' @param object An object of class `fastDR`.
#' @param ... Additional arguments. Currently ignored.
#'
#' @details
#' `summary.fastDR()` prints the estimated treatment effects and the balance
#' table stored on the fitted object. If the object was created with
#' `ps.only = TRUE`, treatment effects are not available and the method reports
#' that no treatment effects were estimated.
#'
#' The balance table compares control and treatment distributions for the
#' covariates used in the propensity score model. See [fastDR()] for details on
#' the returned object components.
#'
#' @return `summary.fastDR()` is called for its side effect of printing model
#'   results and balance diagnostics. It returns `NULL` invisibly.
#'
#' @author Greg Ridgeway \email{gridge@upenn.edu}
#' @seealso [fastDR()], [print.fastDR()]
#'
#' @keywords models
#' @method summary fastDR
#' @export
summary.fastDR <- function(object, ... )
{
   if(!inherits(object, "fastDR"))
   {
      stop("object must be a fastDR object, typically produced from a call to fastDR")
   }

   if(is.null(object$effects))
   {
      cat("fastDR object has no estimated treatment effects. Most likely the call to fastDR had ps.only=TRUE")
   } else
   {
      cat("Results\n")
      cat("Estimand: ", object$estimand, "\n")
      print(object$effects)
   }

   if(is.null(object$balance.tab))
   {
      stop("fastDR object is missing the balance table. Perhaps the model did not run properly.")
   } else
   {
      cat("Balance table\n")
      print(object$balance.tab)
   }
}

#' Weighted Kolmogorov-Smirnov Statistic
#'
#' Computes a weighted Kolmogorov-Smirnov statistic to measure the difference in
#' the marginal feature distributions between the treatment and control cases.
#'
#' @param x A vector of numeric measurements.
#' @param z A vector of 0/1 indicators indicating group membership.
#' @param w A vector of weights.
#'
#' @details
#' This function is used in the propensity score model building step to assess
#' the balance between the treatment and control cases.
#'
#' @return The Kolmogorov-Smirnov statistic, the largest difference between
#'   treatment and control groups' empirical cumulative distribution functions.
#'
#' @references
#' A. Kolmogorov (1933). "Sulla determinazione empirica di una legge di
#' distribuzione," \emph{Giornale dell'Istituto Italiano degli Attuari}
#' 4:83-91.
#'
#' N. Smirnov (1948). "Table for estimating the goodness of fit of empirical
#' distributions," \emph{Annals of Mathematical Statistics} 19(2):279-281.
#'
#' @author Greg Ridgeway \email{gridge@upenn.edu}
#' @seealso [fastDR()]
#'
#' @examples
#' y     <- c(rnorm(100, 0, 1), rnorm(100, 0.5, 1))
#' treat <- rep(0:1, each = 100)
#' w     <- 1 / c(pnorm(y[1:100], 0, 1), pnorm(y[101:200], 0.5, 1))
#' ks(x = y, z = treat, w = w)
#'
#' @keywords univar
#' @export
ks <- function(x,z,w)
{
    w[z==1] <-  w[z==1]/sum(w[z==1])
    w[z==0] <- -w[z==0]/sum(w[z==0])
    ind <- order(x)
    cumv <- abs(cumsum(w[ind]))
    cumv <- cumv[diff(x[ind]) != 0]
    return(ifelse(length(cumv) > 0, max(cumv), 0))
}

.fastDR_paired_prediction_cov_sum <- function(glm1, data, t.var)
{
   data0 <- data
   data1 <- data
   data0[, t.var] <- 0
   data1[, t.var] <- 1
   X0 <- model.matrix(delete.response(terms(glm1)),
                      data=data0,
                      contrasts.arg=glm1$contrasts,
                      xlev=glm1$xlevels)
   X1 <- model.matrix(delete.response(terms(glm1)),
                      data=data1,
                      contrasts.arg=glm1$contrasts,
                      xlev=glm1$xlevels)
   X0 <- X0[, names(coef(glm1)), drop=FALSE]
   X1 <- X1[, names(coef(glm1)), drop=FALSE]
   eta0 <- as.numeric(X0 %*% coef(glm1))
   eta1 <- as.numeric(X1 %*% coef(glm1))
   dmu0 <- glm1$family$mu.eta(eta0)
   dmu1 <- glm1$family$mu.eta(eta1)
   sum(rowSums((X0 %*% vcov(glm1)) * X1) * dmu0 * dmu1)
}

#' Fast Doubly Robust Estimation
#'
#' Doubly robust treatment effect estimation with non-parametric propensity
#' score estimates.
#'
#' @param form.list A list of formulas giving the components of the model,
#'   outcomes, treatment indicator, covariates, observation weights, and
#'   observation keys. Consider using [make.fastDR.formula()] to assist in
#'   creating `form.list`.
#' @param data A data frame containing the variables in the model.
#' @param y.dist The distribution assumed for the outcome. These should be
#'   specified as character strings, for example `"gaussian"`,
#'   `"quasibinomial"`, or `"quasipoisson"`. Using `"binomial"` will cause
#'   inconsequential warning messages since the model uses non-integer
#'   propensity score weights. `"quasibinomial"` avoids this. If `y.form` has
#'   multiple outcomes then `y.dist` should be a vector of strings the same
#'   length as the number of outcomes in `y.dist`.
#' @param estimand Either `"ATT"` (default) for the average treatment effect on
#'   the treated or `"ATE"` for the average treatment effect.
#' @param n.trees The total number of trees to fit for the generalized boosted
#'   model. This is equivalent to the number of iterations and the number of
#'   basis functions in the additive expansion. The default is most likely
#'   appropriate.
#' @param interaction.depth The maximum depth of variable interactions. 1
#'   implies an additive model, 2 implies a model with up to 2-way interactions,
#'   etc.
#' @param shrinkage A shrinkage parameter applied to each tree in the
#'   generalized boosted model expansion. Also known as the learning rate or
#'   step-size reduction. The default is most likely appropriate.
#' @param verbose If `TRUE`, `fastDR()` will print out progress and performance
#'   indicators.
#' @param ps.only If `TRUE`, `fastDR()` will skip the DR step and just return
#'   the propensity scores. This is useful if you need to do non-standard
#'   changes to the regression model, like eliminating some covariates.
#' @param keepGLM If `TRUE`, keep a copy of the `glm` objects produced in the
#'   outcome model step. These models can become large with multiple outcomes
#'   and large datasets.
#' @param smooth.lm A numeric value for penalizing the size of the covariates in
#'   the DR step.
#' @param par_details A list containing the number of threads and
#'   `array_chunk_size` to be passed to `gbmt` for running several gbm
#'   computations in parallel, for example gradient calculation and split
#'   selection. Most easily set using `gbmParallel`.
#'
#' @details
#' The `form.list` has the following components.
#'
#' `y.form` is a formula specifying the outcomes, for example `~y`. Multiple
#' outcomes may be listed, for example `~y1 + y2`; `fastDR()` will conduct a
#' separate analysis for each of them.
#'
#' `t.form` is a formula specifying the treatment indicator, for example
#' `~treat`. The formula can include only one binary treatment coded as 0/1.
#'
#' `x.form` is a formula specifying the potential confounding variables that
#' will be included in the propensity score model and in the outcome model, for
#' example `~X1 + X2 + X2`.
#'
#' `weights.form` is an optional formula specifying observation weights, such as
#' sampling weights. Sampling weights cannot be missing, negative, or equal to
#' zero.
#'
#' `key.form` is a formula giving observation IDs, for example `~caseID`. This
#' is required in order to make sure that after returning results, the propensity
#' score weights can be correctly aligned with cases by the user.
#'
#' `fastDR()` estimates either the average treatment effect on the treated (ATT)
#' or the average treatment effect (ATE).
#'
#' `fastDR()` uses a generalized boosted model to estimate a propensity score
#' model (McCaffrey, Ridgeway, and Morral, 2004). It uses those propensity
#' scores as weights in a weighted generalized linear model to produce a doubly
#' robust estimate of the treatment effect.
#'
#' To reuse the propensity score weights in other analyses, make sure that the
#' code aligns the right weight with the right case. Use the key in the data and
#' the key attached to the weights. For example,
#' `mydata$psw <- myfastDR$w[match(mydata$observationID, names(myfastDR$w))]`.
#'
#' For ATT, treated cases receive weight equal to 1 and comparison cases receive
#' weight p/(1-p). Control cases with factor levels that do not appear among
#' treated cases are dropped from the analysis. For ATE, treated cases receive
#' weight equal to 1/p and comparison cases receive weight 1/(1-p). Cases in
#' either group with factor levels that do not appear in the opposite group are
#' dropped from the analysis, and `fastDR()` issues a warning when this occurs.
#'
#' The effective sample size is the number of independent cases that would yield
#' the same precision as the given weighted sample. It is computed as
#' \deqn{\frac{(\sum_i (1-t_i)w_i)^2}{\sum_i (1-t_i)w_i^2}}
#'
#' When `smooth.lm` is greater than 0, `fastDR()` uses a data augmentation method
#' to impose a penalty on the regression coefficients in the DR step. At a
#' minimum, modest penalties provide numerical stability in the event of
#' correlation among the covariates or certain features being highly predictive
#' of the outcome. `fastDR()` does not penalize the intercept or the treatment
#' indicator. `fastDR()` uses the data augmentation approach described in
#' Greenland and Mansournia (2015). This approach uses ridge regression
#' penalties (N(0,m)) for Gaussian outcomes, log-F(m,m) prior for logistic
#' regression, and gamma(m,m) prior for Poisson regression. The appropriate size
#' of the penalty is subjective, but values between 0.1 and 4.0 are modest
#' penalties for general application.
#'
#' For the DR estimates, `fastDR()` computes exact standard errors when the
#' prediction data used for counterfactual means has at most 5000 rows. This
#' corresponds to at most 2500 treated cases for ATT and at most 2500 analysis
#' cases for ATE. For larger analyses, `fastDR()` takes a random sample of 2500
#' cases for computing the standard errors of the estimates of E(Y0), E(Y1), and
#' E(Y1-Y0). Computing these standard errors requires summarizing a large
#' covariance matrix of predicted values, which can require a lot of memory, but
#' can be efficiently estimated with a sample. Note that there will be some
#' amount of sampling variability; in simulations, the Monte Carlo standard
#' errors are typically off by at most 0.1%.
#'
#' @return `fastDR()` returns a `fastDR` object, a list containing:
#' \item{ks}{The largest KS statistic in the balance table. See `balance.tab`
#'   below. A measure of the quality of the propensity score fit.}
#' \item{ks.un}{The largest KS statistic in the balance table before propensity
#'   score weighting. See `balance.tab.un` below.}
#' \item{balance.tab}{The balance table showing the treatment means and
#'   propensity score weighted control means for each of the terms specified in
#'   `x.form`. For categorical features the table lists each level separately.}
#' \item{balance.tab.un}{The balance table showing the treatment means and
#'   unweighted control means for each of the terms specified in `x.form`.}
#' \item{w}{The propensity score weights used in the analysis. Use `names()` to
#'   extract the key associated with each weight to align with the original
#'   data. All cases in the original data might not necessarily have a
#'   propensity score weight if they were dropped because factor levels were
#'   outside common support. Missing, negative, or zero sampling weights cause
#'   `fastDR()` to stop.}
#' \item{p}{The estimated propensity scores. `p` also maintains the `key` in its
#'   names.}
#' \item{n1, ESS}{The number of cases in the treatment group and the effective
#'   sample size of the control group.}
#' \item{glm.un, glm.ps, glm.dr}{The generalized linear model fits using only
#'   sampling weights (`glm.un`), using propensity score weights (`glm.ps`), and
#'   the model producing doubly robust estimates (`glm.dr`).}
#' \item{effects}{A table showing the results. `effects` includes results for an
#'   unadjusted (`un`), propensity score adjusted (`ps`), and doubly robust
#'   estimate (`dr`). For each, `effects` includes the estimates and standard
#'   errors for E(Y1|t=1), E(Y0|t=1), their difference, and p-value. These
#'   results are all on the scale of the response, not odds ratios or rate
#'   ratios.}
#' \item{y.dist}{Contains the value of `y.dist` used in the call to `fastDR()`.}
#' \item{shrinkage}{The shrinkage parameter used in the propensity score model.}
#'
#' @references
#' D. McCaffrey, G. Ridgeway, and A. Morral (2004). "Propensity score
#' estimation with boosted regression for evaluating adolescent substance abuse
#' treatment," \emph{Psychological Methods} 9(4):403-425.
#'
#' J.H. Friedman (2001). "Greedy Function Approximation: A Gradient Boosting
#' Machine," \emph{Annals of Statistics} 29(5):1189-1232.
#'
#' S. Greenland and M.A. Mansournia (2015). "Penalization, bias reduction, and
#' default priors in logistic and related categorical and survival regressions,"
#' \emph{Statistics in Medicine} 34:3133-3143.
#'
#' S. Greenland, M.A. Mansournia, and D.G. Altman (2014). "Sparse data bias: a
#' problem hiding in plain sight," \emph{The BMJ} 353:i1981.
#'
#' \url{https://gregridgeway.github.io/}
#'
#' @author Greg Ridgeway \email{gridge@upenn.edu}
#' @seealso [gbm3::gbmt()], [gbm3::gbmParallel()]
#'
#' @examples
#' # NHANES example from survey package
#'
#' data(nhanes)
#' nhanes <- nhanes[seq_len(300), ]
#' # add a unique ID to each row
#' nhanes$observationID <- 1:nrow(nhanes)
#' # recode the "treatment" (male) to a 0/1 indicator
#' nhanes$male <- as.numeric(nhanes$RIAGENDR == 1)
#' # make "race" a factor with descriptive labels
#' nhanes$race <- factor(nhanes$race,
#'                       levels = c(1, 2, 3, 4),
#'                       labels = c("Hispanic", "non-Hispanic white",
#'                                  "non-Hispanic black", "other"))
#'
#' dr1 <- fastDR(list(y.form = ~HI_CHOL,          # high cholesterol (over 240mg/dl)
#'                    t.form = ~male,             # compare males to similar females
#'                    x.form = ~race + agecat,    # potential confounders
#'                    weights.form = ~WTMEC2YR,   # sampling weights
#'                    key.form = ~observationID),
#'               data = nhanes,
#'               y.dist = "quasibinomial",   # outcome is 0/1
#'               n.trees = 14,
#'               interaction.depth = 1,
#'               shrinkage = 0.01,
#'               verbose = FALSE,
#'               smooth.lm = 0.1,            # stabilize regression coefs
#'               par_details = gbmParallel(1, 1024)) # just use one core
#'
#' # show balance table
#' round(dr1$balance.tab, 3)
#'
#' # compute percentage difference
#' print(dr1, type = "outcome", model = "dr")
#'
#' @keywords models
#' @export
fastDR <- function(form.list,
                   data,
                   y.dist="gaussian",
                   estimand="ATT",
                   n.trees=3000,
                   interaction.depth=3,
                   shrinkage=0.01,
                   verbose=FALSE,
                   ps.only=FALSE,
                   keepGLM=TRUE,
                   smooth.lm=0,
                   par_details=gbmParallel(1,1024))
{
   y.form       <- form.list$y.form
   if(is.null(y.form))
      stop("form.list is missing y.form")
   t.form       <- form.list$t.form
   if(is.null(t.form))
      stop("form.list is missing t.form")
   x.form       <- form.list$x.form
   if(is.null(x.form))
      stop("form.list is missing x.form")
   weights.form <- form.list$weights.form
   key.form     <- form.list$key.form
   if(is.null(key.form))
      stop("form.list is missing key.form")

   # check that formulas do not use the "."
   if("." %in% all.vars(y.form))
      stop("'.' not allowed in y.form")

   # extract the names of outcome variables
   outcome.y  <- attr(terms(y.form),"term.labels")
   n.outcomes <- length(outcome.y)

   if("." %in% all.vars(x.form))
      warning("You have used a '.' in x.form. This likely included variables that you did not want to be included. If you must use '.' be sure to use '-' to remove the unwanted terms (e.g. ~.-weights-ID-treat-Y)")

   if("." %in% all.vars(t.form))
      stop("'.' not allowed in t.form")
   if(length(all.vars(t.form))>1)
      stop("t.form should be a formula with only one treatment indicator (ex. t.form=~treat)")

   if("." %in% all.vars(key.form))
      stop("'.' not allowed in key.form")
   if(length(all.vars(key.form))>1)
      stop("key.form should have only one term")

   if("." %in% all.vars(weights.form))
      stop("'.' not allowed in weights.form")
   if(length(all.vars(weights.form))>1)
      stop("weights.form should have only one term")

   if(any(!is.character(y.dist)))
      stop("Outcome distribution, y.dist, should be given as a character string (ex. \"gaussian\", with the quotes)")
   if(length(y.dist)==1) y.dist <- rep(y.dist,n.outcomes)
   if(length(y.dist)!=n.outcomes)
      stop("Length of y.dist must be 1 or be the same as the number of outcomes")
   if(any(y.dist=="binomial"))
   {
      y.dist[y.dist=="binomial"] <- "quasibinomial"
      warning("Using 'quasibinomial' instead of 'binomial'")
   }
   if(any(y.dist=="poisson"))
   {
      y.dist[y.dist=="poisson"] <- "quasipoisson"
      warning("Using 'quasipoisson' instead of 'poisson'")
   }
   
   if(!(estimand %in% c("ATT","ATE")))
   {
      stop("estimand must either be 'ATT' or 'ATE'")
   }
   
   if(n.trees<14)
   {
      n.trees <- 14
      warning("n.trees set less than 14. fastDR reset n.trees to 14.")
   }

   # need to use a variable called w and samp.w
   if("w" %in% names(data))
      stop("dataset cannot have a variable names 'w'. fast.dr() needs to use 'w' to store weights")
   if("samp.w" %in% names(data))
      stop("dataset cannot have a variable names 'samp.w'. fast.dr() needs to use 'samp.w' to store sampling weights")

   # check y.form, t.form, and x.form mutually exclusive
   a <- intersect(attr(terms(y.form),"term.labels"),
                  attr(terms(t.form),"term.labels"))
   if(length(a)>0)
      stop("Treatment indicator in t.form is also included in y.form: ",a)
   a <- intersect(attr(terms(t.form),"term.labels"),
                  attr(terms(x.form,data=data),"term.labels"))
   if(length(a)>0)
      stop("Treatment indicator in t.form should not be included in x.form: ",a)
   a <- intersect(attr(terms(y.form),"term.labels"),
                  attr(terms(x.form,data=data),"term.labels"))
   if(length(a)>0)
      stop("Outcome in y.form should not be included in x.form: ",a)

   # extract observation key (case ID)
   if(is(key.form)[1]=="formula")
   {
      key <- model.frame(key.form,data, na.action=na.pass)[,1]
   }
   if(any(is.na(key)))
      stop("missing values in key. key cannot have missing values")
   if(length(key)!=nrow(data))
      stop("The length of key must be the same as the number of rows in data. key should be a unique case ID")
   if(length(key)!=length(unique(key)))
      stop("key has duplicates. key should be a unique case ID")

   # extract dataset
   data0 <- cbind(model.frame(y.form, data, na.action=na.pass),
                  model.frame(t.form, data, na.action=na.pass),
                  model.frame(x.form, data, na.action=na.pass))
   
   # get treatment indicator and check that they are all 0/1
   i.treat <- model.frame(t.form,data, na.action=na.pass)[,1]
   if(!all(i.treat %in% 0:1))
      stop("The treatment indicator specified in t.form does not only take values 0 or 1. The treatment indicator must be a 0 or a 1")
   # expressions that can be eval'd in subset() to select treat/control cases
   subsetExpr0 <- parse(text=paste(as.character(t.form[[2]]),"==0"))
   subsetExpr1 <- parse(text=paste(as.character(t.form[[2]]),"==1"))

   # extract the variables names in X
   match.vars <- names(data0)[-(1:(n.outcomes+1))]
   match.vars.NA.index <- NULL

   # extract observation weights
   if(is.null(weights.form))
   {
      data0$samp.w <- rep(1,nrow(data))
   } else
   {
      if(is(weights.form)[1]!="formula")
         stop("weights parameter should be a formula of the form ~weight, where 'weight' is the variable in data containing the observation weights")
      warning("SEs for DR estimates with sample weighted data might not completely account for sampling weights")
      data0$samp.w <- model.frame(weights.form,data, na.action=na.pass)[,1]
   }
   if(any(is.na(data0$samp.w)))
      stop("missing values in weights. weights cannot have missing values")
   if(any(data0$samp.w<0))
      stop("Some observation weights are negative")
   if(any(data0$samp.w==0))
      stop("Some observation weights are 0. Drop cases with sampling weights equal to 0 before calling fastDR")

   # create missing indicators, median impute numeric features
   for(xj in match.vars)
   {
      i <- which(is.na(data0[,xj]))
      if(!is.factor(data0[[xj]]) && (length(i)>10))
      {  # create missing indicators if there are at least 10 NAs
         #   check not completely correlated with other missing indicators
         a <- rep(0, nrow(data0))
         a[i] <- 1
         corNA <- cor(a,data0[,match.vars.NA.index])
         if(all(corNA!=1))
         {
            data0 <- cbind(data0, a)
            names(data0)[ncol(data0)] <- paste(xj,".NA",sep="")
            match.vars.NA.index <- c(match.vars.NA.index, ncol(data0))
         }
      }
      if(length(i)>0)
      {
         if(is.factor(data0[,xj]))
         {  # set NA to be a valid factor level
            data0[[xj]] <- factor(data0[[xj]],exclude=NULL)
            levels(data0[[xj]])[is.na(levels(data0[[xj]]))] <- "<NA>"
         }
         else # numeric
         {  # impute the median
            data0[i,xj] <- median(data0[,xj],na.rm=TRUE)
         }
      }
   }

   # in case the new missing indicators duplicate existing variable names
   names(data0) <- make.unique(names(data0))
   match.vars.NA <- names(data0)[match.vars.NA.index]

   # include missing indicators in x.form
   if(length(match.vars.NA)>0)
   {
      x.form <- update(x.form, formula(paste("~.+",
                       paste(match.vars.NA,collapse="+"),sep="")),
                       data=data)
      match.vars <- c(match.vars,match.vars.NA)
   }

   ps.form <- formula(paste(as.character(t.form[[2]]),"~",
                            as.character(x.form[2])))

   # drop cases with factor levels outside common support
   k <- rep(TRUE,nrow(data0))
   for(xj in match.vars[sapply(data0[1,match.vars],is.factor)])
   {
      treat.levels <- unique(subset(data0, eval(subsetExpr1))[[xj]])
      cntrl.levels <- unique(subset(data0, eval(subsetExpr0))[[xj]])
      k.cntrl <- is.na(data0[,xj]) | (data0[,xj] %in% treat.levels)
      k.treat <- is.na(data0[,xj]) | (data0[,xj] %in% cntrl.levels)
      if(estimand=="ATE")
      {
         dropped.cntrl <- sum(k & with(data0, eval(subsetExpr0)) & !k.cntrl)
         dropped.treat <- sum(k & with(data0, eval(subsetExpr1)) & !k.treat)
         if((dropped.cntrl+dropped.treat)>0)
         {
            warning("For ATE, dropping ", dropped.cntrl,
                    " control and ", dropped.treat,
                    " treatment cases with factor levels in ", xj,
                    " that do not appear in the opposite group")
         }
         k <- k & k.cntrl & k.treat
      } else
      {
         k <- k & k.cntrl
      }
   }

   data0   <- data0[k,]
   key     <- key[k]

   # drop empty levels
   for(j in match.vars)
   {
      if(is.factor(data0[,j]))
      {
         data0[,j] <- factor(data0[,j])
         # eliminate one-level factor variables from DR, svyglm would fail
         if(nlevels(data0[,j])<=1)
         {
            warning("Dropping factor variable with one-level: ",j)
            x.form <- update(x.form,formula(paste0("~ . -",j)))
         }
      }
   }
   results <- list()

### BEGIN PROPENSITY SCORE ESTIMATION ###
   if(verbose)
      cat("Fitting generalized boosted model for propensity scores\n")
   converged <- FALSE
   while(!converged)
   {
      if(verbose) cat("shrinkage:",round(shrinkage,4),"\n")
      gbm1 <- gbmt(ps.form,
                   distribution=gbm_dist("Bernoulli"),
                   data=data0,
                   weights=data0$samp.w,
                   train_params=
                      training_params(num_trees=n.trees,
                                      interaction_depth=interaction.depth,
                                      shrinkage=shrinkage,
                                      bag_fraction=1.0,
                                      num_train=nrow(data0),
                                      num_features=length(match.vars)),
                   is_verbose=verbose,
                   keep_gbm_data=FALSE,
                   par_details=par_details)

      iters <- unique(c(0, as.integer(round(seq(trunc(n.trees*0.3),
                                                n.trees,
                                                length.out=11)))))
      p <- predict(gbm1,newdata=data0,n.trees=iters,type="response")

      if(verbose)
         cat("Assessing balance...\n")
      best.bal <- 1
      for(j in 1:ncol(p))
      {
         i.cntrl <- with(data0, eval(subsetExpr0))
         i.treat <- with(data0, eval(subsetExpr1))
         data0$w <- 1
         if(estimand=="ATT")
         {
            data0$w[i.cntrl] <- p[i.cntrl,j]/(1-p[i.cntrl,j])
         } else # estimand=="ATE"
         {
            data0$w[i.treat] <- 1/p[i.treat,j]
            data0$w[i.cntrl] <- 1/(1-p[i.cntrl,j])
         }
         
         # weights should be sampling weight*PSW
         # G. Ridgeway, S. Kovalchik, B.A. Griffin, and M.U. Kabeto (2015). 
         #   "Propensity score analysis with survey weighted data," Journal of
         #   Causal Inference 3(2):237-249
         data0$w <- with(data0, samp.w*w) 

         # normalize the weights to be max 1.0 within treatment group
         data0$w[i.treat] <- data0$w[i.treat]/max(data0$w[i.treat])
         data0$w[i.cntrl] <- data0$w[i.cntrl]/max(data0$w[i.cntrl])

         bal <- NULL
         for(x in match.vars)
         {
            if(is.factor(data0[,x]))
            {
               a <- cbind(sapply(split(data0$w[i.cntrl], data0[i.cntrl,x]), sum),
                          sapply(split(data0$w[i.treat], data0[i.treat,x]), sum))
               a <- t(t(a)/colSums(a))
               a <- cbind(a, a[,2]-a[,1])
               rownames(a) <- paste(x,":",rownames(a),sep="")
            } else
            {
               a <- c(weighted.mean(data0[i.cntrl,x], data0$w[i.cntrl], na.rm=TRUE),
                      weighted.mean(data0[i.treat,x], data0$w[i.treat] ,na.rm=TRUE))
               a <- c(a, ks(data0[,x],i.treat,data0$w))
               a <- matrix(a,nrow=1)
               rownames(a) <- x
            }
            bal <- rbind(bal,a)
         }

         this.bal <- max(abs(bal[,3]))
         if(is.na(this.bal))
         {
            this.bal <- Inf
         }
         if((best.bal>0.002) &&    # don't refine further if best KS<0.002
            (this.bal < best.bal))
         {
            best.bal <- this.bal
            results$ks <- best.bal
            results$best.iter <- iters[j]
            results$balance.tab <- bal
            results$w <- data0$w
            results$p <- p[,j]
         }
         if(j==1) # unweighted analysis
         {
            results$ks.un <- this.bal
            results$balance.tab.un <- bal
         }
      }

      # add key to the names of the weights so user can match weights to cases
      names(results$w) <- key
      names(results$p) <- key
      results$shrinkage <- shrinkage
      if(verbose)
      {
         cat("Best number of iterations:",results$best.iter,"\n")
         cat("Max KS:",results$ks,"\n")
      }
      converged <-
         with(results, (ks<0.005) ||
                       (best.iter>n.trees*0.5 && best.iter<n.trees*0.9))
      if(verbose && (results$ks<0.005))
         cat("Propensity score model converged with KS<0.005\n")
      if(results$best.iter==0)
      {
         converged <- TRUE
         warning("Best iteration is 0 in propensity score model. Might signal a problem with the data")
      }
      if(shrinkage*results$best.iter/(n.trees*0.75) > 1)
      {
         converged <- TRUE
         warning("Projected shrinkage would exceed 1.0 if we keep going. This is probably the best propensity score model gbm can produce.")
      } else
      {
         shrinkage <- shrinkage*results$best.iter/(n.trees*0.75)
      }
   }

   i.cntrl <- with(data0, which(eval(subsetExpr0)))
   i.treat <- with(data0, which(eval(subsetExpr1)))
   results$n1  <- length(i.treat)
   w0 <- results$w[i.cntrl]
   results$ESS <- sum(w0)^2/sum(w0^2)
   colnames(results$balance.tab.un) <- c("control","treatment","KS")
   colnames(results$balance.tab)    <- c("control","treatment","KS")
   if(verbose)
      cat("Propensity score estimation complete\n")
### END PROPENSITY SCORE ESTIMATION ###

   if(ps.only)
   {
      class(results) <- "fastDR"
      return(results)
   }

### BEGIN OUTCOME ANALYSIS ###
   if(keepGLM)
   {
      results$glm.un <- vector("list",length(outcome.y))
      results$glm.ps <- vector("list",length(outcome.y))
      results$glm.dr <- vector("list",length(outcome.y))
   }
   results$z      <- rep(NA,length(outcome.y))

   # make sure data0 has the best prop score weights
   data0$w <- results$w
   data0$samp.w[i.cntrl] <- data0$samp.w[i.cntrl]/max(data0$samp.w[i.cntrl])
   data0$samp.w[i.treat] <- data0$samp.w[i.treat]/max(data0$samp.w[i.treat])
   sdesign.un <- svydesign(ids=~1, weights=~samp.w, data=data0)
   sdesign.w  <- svydesign(ids=~1, weights=~w,      data=data0)

   # for storing the results
   results$effects <- vector("list",length(outcome.y))
   names(results$effects) <- attr(terms(y.form),"term.labels")
   
   if(verbose) cat("Fitting outcome regression models...")
   for(i.y in 1:length(outcome.y))
   {
      results$effects[[i.y]] <-
         data.frame(E.y1 =rep(0,3), E.y0 =rep(0,3),
                    se.y1=rep(0,3), se.y0=rep(0,3),
                    TE   =rep(0,3), se.TE=rep(0,3),
                    p    =rep(0,3))
      rownames(results$effects[[i.y]]) <- c("un","ps","dr")
      
      ps.form <- formula(paste(outcome.y[i.y],"~",
                               as.character(t.form[[2]])))
      if(verbose) cat(outcome.y[i.y],"...")

      ### unweighted analysis ###
      # must use "substitute" to inject arguments into svyglm, scoping issue
      glm1 <- substitute(svyglm(formula = ps.form, 
                                design  = sdesign.un,
                                family  = y.dist[i.y],
                                rescale = FALSE))
      glm1 <- eval(glm1)
      if(keepGLM) results$glm.un[[i.y]] <- glm1
      
      a <- svymean(reformulate(attr(terms(y.form), "term.labels")[i.y]),
                   design=subset(sdesign.un, eval(subsetExpr1)),
                   na.rm=TRUE)
      results$effects[[i.y]]$E.y1[]   <- as.numeric(a)
      results$effects[[i.y]]$se.y1[]  <- sqrt(vcov(a))
      
      a <- svymean(reformulate(attr(terms(y.form), "term.labels")[i.y]),
                   design=subset(sdesign.un, eval(subsetExpr0)),
                   na.rm=TRUE)
      results$effects[[i.y]]$E.y0[1]  <- as.numeric(a)
      results$effects[[i.y]]$se.y0[1] <- sqrt(vcov(a))

      # grab a generic treatment row and control row
      dataTwoRows01 <- rbind(subset(data0, eval(subsetExpr0))[1,],
                             subset(data0, eval(subsetExpr1))[1,])
      y.hat0 <- predict(glm1, newdata=dataTwoRows01, type="response", vcov=TRUE)
      
      results$effects[[i.y]]$TE[1]    <- as.numeric(t(c(-1,1)) %*% y.hat0)
      results$effects[[i.y]]$se.TE[1] <- as.numeric(sqrt(t(c(-1,1)) %*% vcov(y.hat0) %*% c(-1,1)))
      results$effects[[i.y]]$p[1]     <- coef(summary(glm1))[2,4]
       
      ### propensity score ###
      glm1 <- substitute(svyglm(formula = ps.form,
                                design  = sdesign.w,
                                family  = y.dist[i.y],
                                rescale = FALSE))
      glm1 <- eval(glm1)
      if(keepGLM) results$glm.ps[[i.y]] <- glm1

      a <- svymean(reformulate(attr(terms(y.form), "term.labels")[i.y]),
                   design=subset(sdesign.w, eval(subsetExpr0)),
                   na.rm=TRUE)
      results$effects[[i.y]]$E.y0[2]  <- as.numeric(a)
      results$effects[[i.y]]$se.y0[2] <- sqrt(vcov(a))
      if(estimand=="ATE")
      {
         a <- svymean(reformulate(attr(terms(y.form), "term.labels")[i.y]),
                      design=subset(sdesign.w, eval(subsetExpr1)),
                      na.rm=TRUE)
         results$effects[[i.y]]$E.y1[2]  <- as.numeric(a)
         results$effects[[i.y]]$se.y1[2] <- sqrt(vcov(a))
      }
      
      y.hat0 <- predict(glm1, newdata=dataTwoRows01, type="response", vcov=TRUE)
      results$effects[[i.y]]$TE[2]    <- t(c(-1,1)) %*% y.hat0
      results$effects[[i.y]]$se.TE[2] <- sqrt(t(c(-1,1)) %*% vcov(y.hat0) %*% c(-1,1))
      results$effects[[i.y]]$p[2]     <- coef(summary(glm1))[2,4]
      
      ### DR ###
      nMissingOutcome <- sum(is.na(data0[[outcome.y[i.y]]]))
      if(nMissingOutcome > 0)
         warning("Outcome ",outcome.y[i.y]," has ",nMissingOutcome," missing values. They are included in the propensity score stage but dropped from the regression step")
      
      #    create intercept here... need control of it for prior penalty
      data.mx <- model.matrix(formula(paste("~w+",outcome.y[i.y],"+",
                                            as.character(t.form[[2]]),"+",
                                            as.character(x.form[2]))),
                              data=data0)
      colnames(data.mx)[1] <- "Intercept"
      
      # put penalty on regression terms
      # S. Greenland M.A.  Mansournia (2015). "Penalization, bias reduction, and 
      # default priors in logistic and related categorical and survival 
      # regressions. Statistics in Medicine 34:3133-3143. 10.1002/sim.6537.
      
      # S. Greenland, M.A. Mansournia, D.G. Altman (2016). "Sparse data bias: a 
      # problem hiding in plain sight," BMJ 352:i1981 10.1136/bmj.i1981.
      if(smooth.lm>0)
      {
         a <- cbind(0, smooth.lm, 0, 0, diag(1, nrow = ncol(data.mx)-4))
         data.mx <- rbind(data.mx, a)
         if(y.dist[[i.y]]=="quasibinomial") # for log-F distribution need two rows
         {
            a[,3] <- 1 # set outcome to 1
            data.mx <- rbind(data.mx, a)
         }
      }
      
      data.mx <- data.frame(data.mx, row.names = 1:nrow(data.mx))
      sdesign.wexp <- svydesign(ids=~1, weights=~w, data=data.mx)
      # Static check hint; subset() evaluates Intercept from data.mx.
      Intercept <- NULL
      
      dr.form <- formula(paste(outcome.y[i.y],"~-1+Intercept+",
                               paste(names(data.mx)[-(1:3)],collapse="+")))
      converged <- FALSE
      while(!converged)
      {
         glm1 <- substitute(svyglm(formula = dr.form,
                                   design  = sdesign.wexp,
                                   family  = y.dist[i.y],
                                   rescale = FALSE))
         glm1 <- eval(glm1)
         converged <- all(!is.na(glm1$coefficients))
         if(!converged)
         {
            to.drop <- names(glm1$coefficients[is.na(glm1$coefficients)])
            to.drop <- setdiff(to.drop,as.character(t.form[[2]]))
            if(all(to.drop==as.character(t.form[[2]])))
            {
               glm1 <- NA
               converged <- TRUE
            } else
            {
               cat("Dropping",to.drop,"\n")
               dr.form <- update(dr.form,formula(paste(".~.",paste0("-",to.drop,collapse=""))))
            }
         }
      }
      
      if(keepGLM) results$glm.dr[[i.y]] <- glm1
      results$z[i.y] <- coef(summary(glm1))[2,"t value"]
      if(FALSE) # transform from t to z for better FDR calculation
      {
       results$z[i.y] <- qnorm(pt(coef(summary(glm1))[2,"t value"],
                                  glm1$df.residual))
      }
      
      # collect DR statistics
      results$effects[[i.y]]$p[3] <- coef(summary(glm1))[2,4]
      
      # only real cases, not those for shrinking beta
      if(estimand=="ATT")
      {
         pred.data <- subset(data.mx, Intercept==1 & eval(subsetExpr1))
      } else # estimand=="ATE"
      {
         pred.data <- subset(data.mx, Intercept==1)
      }
      n <- nrow(pred.data)
      t.var <- as.character(t.form[2])
      a <- rbind(pred.data,pred.data)
      a[1:n, t.var] <- 0 # recode cases as control
      a[(n+1):(2*n), t.var] <- 1 # recode cases as treated
      y.hat0 <- predict(glm1,
                        newdata=a,
                        type="response")
      results$effects[[i.y]]$E.y0[3] <- mean(y.hat0[1:n])
      if(estimand=="ATE")
      {
         results$effects[[i.y]]$E.y1[3] <- mean(y.hat0[(n+1):(2*n)])
      }
      results$effects[[i.y]]$TE[3]   <- with(results$effects[[i.y]], E.y1[3]-E.y0[3])
      
      if(sign(coef(glm1)[2]) != 
         sign(results$effects[[i.y]]$TE[3]))
      {
         warning("Outcome regression model treatment coefficient has a different sign than the estimated treatment effect. That is a little unusual and might need a closer look.")   
      }
      
      if(nrow(a)<=5000)
      {
         if(verbose) cat("Exact SE calculation\n")
         u <- cbind(rep(1:0,each=n),         # for SE(EY0)
                    rep(0:1,each=n),         # for SE(EY1)
                    rep(c(-1,1),each=n))/n   # for SE(EY1-EY0)
         y.hat0 <- predict(glm1, 
                           newdata=a, 
                           type="response", 
                           vcov=TRUE)
         se.TE <- sqrt(diag(t(u) %*% vcov(y.hat0) %*% u))
      } else
      {
         n0 <- 2500
         Vdiag <- as.numeric(vcov(y.hat0)) # vcov will return a vector here
         VdiagE0 <- sum(Vdiag[1:n])
         VdiagE1 <- sum(Vdiag[(n+1):(2*n)])
         VdiagY0Y1 <- .fastDR_paired_prediction_cov_sum(glm1,
                                                        pred.data,
                                                        t.var)
         VdiagTE <- VdiagE0 + VdiagE1 - 2*VdiagY0Y1
         
         a <- pred.data[sample(1:n, size=n0),]
         a <- rbind(a,a)
         a[1:n0, t.var] <- 0 # recode sampled cases as control
         a[(n0+1):(2*n0), t.var] <- 1 # recode sampled cases as treated
         y.hat0 <- predict(glm1, 
                           newdata=a, 
                           type="response", 
                           vcov=TRUE)
         
         # se(EY0)
         V <- vcov(y.hat0)[1:n0,1:n0]
         VoffdiagE0 <- sum(V)-sum(diag(V))
         VoffdiagE0 <- n*(n-1)*VoffdiagE0/(n0*(n0-1))
         VEY0 <- VdiagE0 + VoffdiagE0
         # se(EY1)
         V <- vcov(y.hat0)[(n0+1):(2*n0), (n0+1):(2*n0)]
         VoffdiagE1 <- sum(V)-sum(diag(V))
         VoffdiagE1 <- n*(n-1)*VoffdiagE1/(n0*(n0-1))
         VEY1 <- VdiagE1 + VoffdiagE1
         
         # se(EY1-EY0)
         V <- vcov(y.hat0)
         Vcross <- V[1:n0, (n0+1):(2*n0), drop=FALSE]
         VoffdiagY0Y1 <- sum(Vcross) - sum(diag(Vcross))
         VoffdiagY0Y1 <- n*(n-1)*VoffdiagY0Y1/(n0*(n0-1))
         VTE <- VdiagTE + VoffdiagE0 + VoffdiagE1 - 2*VoffdiagY0Y1
         se.TE <- sqrt(c(VEY0, VEY1, VTE))/n
      }
      results$effects[[i.y]]$se.y0[3] <- se.TE[1]
      results$effects[[i.y]]$se.y1[3] <- se.TE[2]
      results$effects[[i.y]]$se.TE[3] <- se.TE[3]
   }
   if(verbose) cat("\nOutcome regression models complete\n")
   ### END OUTCOME ANALYSIS ###

   results$y.dist <- y.dist
   results$estimand <- estimand

   class(results) <- "fastDR"
   return(results)
}
