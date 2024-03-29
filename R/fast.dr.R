# Notes:
#   - fix treatment cases weighted with sampling weights


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
      if(!(weights.var %in% names(data)))
         warning("weights.var does not appear in data")
      if(!(key.var %in% names(data)))
         warning("key.var does not appear in data")
      if(!(key.var %in% names(data)))
         warning("key.var does not appear in data")
      if(x.vars==".")
         x.vars <- setdiff(names(data),c(y.vars,t.var,weights.var,key.var))
   }
   if(is.null(data) && (x.vars="."))
      stop("If x.vars given as '.' then data must be given")

   return(list(
      y.form      =paste("~",paste(y.vars,collapse="+"),sep=""),
      t.form      =paste("~",t.var,sep=""),
      x.form      =paste("~",paste(x.vars,collapse="+"),sep=""),
      weights.form=paste("~",weights.var,sep=""),
      key.form    =paste("~",key.var,sep="")))
}

print.fastDR <- function(x, type="outcome", model="dr",... )
{
   if(class(x)!="fastDR")
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

summary.fastDR <- function(object, ... )
{
   if(class(object)!="fastDR")
   {
      stop("object must be a fastDR object, typically produced from a call to fastDR")
   }

   if(is.null(object$effects))
   {
      cat("fastDR object has no estimated treatment effects. Most likely the call to fastDR had ps.only=TRUE")
   } else
   {
      cat("Results\n")
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

ks <- function(x,z,w)
{
    w[z==1] <-  w[z==1]/sum(w[z==1])
    w[z==0] <- -w[z==0]/sum(w[z==0])
    ind <- order(x)
    cumv <- abs(cumsum(w[ind]))
    cumv <- cumv[diff(x[ind]) != 0]
    return(ifelse(length(cumv) > 0, max(cumv), 0))
}

fastDR <- function(form.list,
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
      key <- model.frame(key.form,data)[,1]
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
      data0$samp.w <- model.frame(weights.form,data)[,1]
   }
   if(any(is.na(data0$samp.w)))
      stop("missing values in weights. weights cannot have missing values")
   if(any(data0$samp.w<0))
      stop("Some observation weights are negative")

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

   # drop control cases with factor levels that treatment cases do not have
   #   saves computation time if eliminate control cases with prop score=0
   k <- rep(TRUE,nrow(data0))
   for(xj in match.vars[sapply(data0[1,match.vars],is.factor)])
   {
      a <- subset(data0, eval(subsetExpr1))[[xj]]
      k <- k & (is.na(data0[,xj]) | (data0[,xj] %in% a))
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

      iters <- c(0,seq(trunc(n.trees*0.3),n.trees,length=11))
      p <- predict(gbm1,newdata=data0,n.trees=iters,type="response")

      if(verbose)
         cat("Assessing balance...\n")
      best.bal <- 1
      for(j in 1:ncol(p))
      {
         data0$w <- 1
         i.cntrl <- with(data0, eval(subsetExpr0))
         i.treat <- with(data0, eval(subsetExpr1))
         data0$w[i.cntrl] <- p[i.cntrl,j]/(1-p[i.cntrl,j])
         
         # weights should be sampling weight*PSW
         # G. Ridgeway, S. Kovalchik, B.A. Griffin, and M.U. Kabeto (2015). 
         #   "Propensity score analysis with survey weighted data,” Journal of 
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
   # this produces incorrect estimates with multiple outcomes with missing values
   # means.un1 <- svymean(y.form,
   #                      design=subset(sdesign.un, eval(subsetExpr1)),
   #                      na.rm=TRUE)
   # means.un0 <- svymean(y.form, 
   #                      design=subset(sdesign.un, eval(subsetExpr0)),
   #                      na.rm=TRUE)
   # means.ps0 <- svymean(y.form,
   #                      design=subset(sdesign.w, eval(subsetExpr0)),
   #                      na.rm=TRUE)
   
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
      # regressions. Statistics in Medicine 34:3133–3143. 10.1002/sim.6537.
      
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
      
      # only real treated cases, not those for shrinking beta
      a <- subset(data.mx, Intercept==1 & eval(subsetExpr1))
      n <- nrow(a)
      a <- rbind(a,a)
      a[1:n, as.character(t.form[2])] <- 0 # recode treated cases as control
      y.hat0 <- predict(glm1,
                        newdata=a,
                        type="response")
      results$effects[[i.y]]$E.y0[3] <- mean(y.hat0[1:n])
      results$effects[[i.y]]$TE[3]   <- with(results$effects[[i.y]], E.y1[3]-E.y0[3])
      
      if(sign(coef(glm1)[2]) != 
         sign(results$effects[[i.y]]$TE[3]))
      {
         warning("Outcome regression model treatment coefficient has a different sign than the estimated treatment effect. That is a little unusual and might need a closer look.")   
      }
      
      if(nrow(a)<=3000)
      {
         u <- cbind(rep(1:0,each=n),         # for SE(EY0)
                    rep(c(-1,1),each=n))/n   # for SE(EY1-EY0)
         y.hat0 <- predict(glm1, 
                           newdata=a, 
                           type="response", 
                           vcov=TRUE)
         se.TE <- sqrt(diag(t(u) %*% vcov(y.hat0) %*% u))
      } else
      {
         VdiagE0 <- sum(vcov(y.hat0)[1:n]) # vcov will return a vector here
         VdiagTE <- sum(vcov(y.hat0))
         
         n0 <- 1500
         a <- subset(data.mx, Intercept==1 & eval(subsetExpr1))
         a <- a[sample(1:nrow(a), size=n0),]
         a <- rbind(a,a)
         a[1:n0, as.character(t.form[2])] <- 0 # recode treated cases as control
         y.hat0 <- predict(glm1, 
                           newdata=a, 
                           type="response", 
                           vcov=TRUE)
         
         # se(EY0)
         V <- vcov(y.hat0)[1:n0,1:n0]
         VoffdiagE0 <- sum(V)-sum(diag(V))
         VoffdiagE0 <- n*(n-1)*VoffdiagE0/(n0*(n0-1))
         VEY0 <- VdiagE0 + VoffdiagE0
         
         # se(ET1-EY0)
         # for the difference variance is
         #   mean(V[1:n,1:n])-2*mean(V[1:n,-(1:n)])+mean(V[-(1:n),-(1:n)])
         V <- vcov(y.hat0)
         VTE <- VdiagTE + n*(n-1)*
            (sum(V[  1:n0,   1:n0]) +
                sum(V[-(1:n0),-(1:n0)]) -
                sum(diag(V)) -
                2*sum(V[1:n0,-(1:n0)])) / (n0*(n0-1))
         se.TE <- sqrt(c(VEY0, VTE))/n
      }
      results$effects[[i.y]]$se.y0[3] <- se.TE[1]
      results$effects[[i.y]]$se.TE[3] <- se.TE[2]
   }
   if(verbose) cat("\nOutcome regression models complete\n")
   ### END OUTCOME ANALYSIS ###

   results$y.dist <- y.dist

   class(results) <- "fastDR"
   return(results)
}
