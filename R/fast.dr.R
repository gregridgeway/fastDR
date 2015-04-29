# Notes:
#   - fix treatment cases weighted with sampling weights


.onAttach <- function(libname, pkgname)
{
   vers <- library(help=fastDR)$info[[1]]
   vers <- vers[grep("Version:",vers)]
   vers <- rev(strsplit(vers," ")[[1]])[1]
   packageStartupMessage(paste("Loaded fastDR",vers))
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
   if(!(type %in% c("outcome","complete")))
      stop("type parameter must be one of 'outcome' (default) or 'complete'")

   if(type=="outcome")
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
         if(x$y.dist[i] %in% c("binomial","quasibinomial"))
         {
            temp <- signif(exp(temp),3)
            cat("odds-ratio (95% CI): ",temp[1],
                " (",temp[2],",",temp[3],")\n",sep="")
         }
         else if(x$y.dist[i] %in% c("poisson","quasipoisson"))
         {
            temp <- signif(exp(temp),3)
            cat("rate-ratio (95% CI): ",temp[1],
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
   cat("Results\n")
   print(object$effects)
   cat("Balance table\n")
   print(object$balance.tab)
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
                   verbose=FALSE)
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
   if(any(y.dist=="binomial"))
      warning("Suggest using 'quasibinomial' instead of 'binomial'")
   if(any(y.dist=="poisson"))
      warning("Suggest using 'quasipoisson' instead of 'poisson'")

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
   i.treat <- i.treat==1

   match.vars <- names(data0)[-(1:2)]
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

   # create missing indicators, and median/mode impute
   for(xj in match.vars)
   {
      i <- which(is.na(data0[,xj]))
      if(length(i)>10)
      {  # create missing indicators if there are at least 10 NAs
         data0 <- cbind(data0, 0)
         data0[i,ncol(data0)] <- 1
         names(data0)[ncol(data0)] <- paste(xj,".NA",sep="")
         match.vars.NA.index <- c(match.vars.NA.index, ncol(data0))
      }
      if(length(i)>0)
      {
         if(is.factor(data0[,xj]))
         {  # impute the mode
            a <- table(data0[,xj])
            data0[i,xj] <- names(a)[which.max(a)]
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
      a <- unique(data0[i.treat,xj])
      k <- k & (is.na(data0[,xj]) | (data0[,xj] %in% a))
   }

   data0   <- data0[k,]
   key     <- key[k]
   i.treat <- i.treat[k]

   # drop empty levels
   for(j in match.vars)
      if(is.factor(data0[,j])) data0[,j] <- factor(data0[,j])

   results <- list()

### BEGIN PROPENSITY SCORE ESTIMATION ###
   if(verbose)
      cat("Fitting generalized boosted model for propensity scores\n")
   converged <- FALSE
   while(!converged)
   {
      if(verbose) cat("shrinkage:",round(shrinkage,4),"\n")
      gbm1 <- gbm(ps.form,
                  data=data0,
                  weights=data0$samp.w,
                  distribution="bernoulli",
                  n.trees=n.trees,
                  interaction.depth=interaction.depth,
                  shrinkage=shrinkage,
                  verbose=verbose,
                  keep.data=FALSE,
                  bag.fraction=1.0)

      iters <- c(0,seq(trunc(n.trees*0.3),n.trees,length=11))
      p <- predict(gbm1,newdata=data0,n.trees=iters,type="response")

      if(verbose)
         cat("Assessing balance...\n")
      best.bal <- 1
      for(j in 1:ncol(p))
      {
         data0$w[i.treat]  <- 1
         data0$w[!i.treat] <- p[!i.treat,j]/(1-p[!i.treat,j])
         data0$w <- with(data0, samp.w*w) # weights should be sampling weight*PSW

         # normalize the weights to have mean 1.0 within treatment
         data0$w[ i.treat] <- data0$w[ i.treat]/mean(data0$w[ i.treat])
         data0$w[!i.treat] <- data0$w[!i.treat]/mean(data0$w[!i.treat])
         
         bal <- NULL
         for(x in match.vars)
         {
            if(is.factor(data0[,x]))
            {
               a <- cbind(sapply(split(data0$w[!i.treat], data0[!i.treat,x]), sum),
                          sapply(split(data0$w[i.treat],  data0[i.treat,x]),  sum))
               a <- t(t(a)/colSums(a))
               a <- cbind(a, a[,2]-a[,1])
               rownames(a) <- paste(x,":",rownames(a),sep="")
            } else
            {
               a <- c(weighted.mean(data0[!i.treat,x],data0$w[!i.treat],na.rm=TRUE),
                      weighted.mean(data0[i.treat,x], data0$w[i.treat] ,na.rm=TRUE))
               a <- c(a, ks(data0[,x],i.treat,data0$w))
               a <- matrix(a,nrow=1)
               rownames(a) <- x
            }
            bal <- rbind(bal,a)
         }

         this.bal <- max(abs(bal[,3]))
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
         cat("Best number of iterations:",results$best.iter,"\n")
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
      shrinkage <- shrinkage*results$best.iter/(n.trees*0.75)
   }

   results$n1  <- sum(i.treat)
   w0 <- results$w[!i.treat]
   results$ESS <- sum(w0)^2/sum(w0^2)
   colnames(results$balance.tab.un) <- c("control","treatment","KS")
   colnames(results$balance.tab)    <- c("control","treatment","KS")
   if(verbose)
      cat("Propensity score estimation complete\n")
### END PROPENSITY SCORE ESTIMATION ###

### BEGIN OUTCOME ANALYSIS ###
   outcome.y <- attr(terms(y.form),"term.labels")
   if(length(y.dist)==1) y.dist <- rep(y.dist,length(outcome.y))
   if(length(y.dist)!=length(outcome.y))
      stop("Length of y.dist must be 1 or be the same as the number of outcomes")
   results$glm.un <- vector("list",length(outcome.y))
   results$glm.ps <- vector("list",length(outcome.y))
   results$glm.dr <- vector("list",length(outcome.y))
   results$z      <- rep(NA,length(outcome.y))

   # make sure data0 has the best prop score weights
   data0$w <- results$w
   sdesign.un <- svydesign(ids=~1,weights=~samp.w,data=data0)
   sdesign.w  <- svydesign(ids=~1,weights=~w,     data=data0)
   if(verbose)
      cat("Fitting outcome regression models...")

   for(i.y in 1:length(outcome.y))
   {
       ps.form <- formula(paste(outcome.y[i.y],"~",
                                as.character(t.form[[2]])))
       if(verbose) cat(outcome.y[i.y],"...")

       # unweighted analysis
       # must use "substitute" to inject arguments into svyglm, scoping issue
       glm1 <- substitute(svyglm(formula=ps.form,design=sdesign.un,
                                 family=y.dist[i.y]))
       glm1 <- eval(glm1)
       results$glm.un[[i.y]] <- glm1

       # propensity score
       glm1 <- substitute(svyglm(formula=ps.form,design=sdesign.w,
                                 family=y.dist[i.y]))
       glm1 <- eval(glm1)
       results$glm.ps[[i.y]] <- glm1

       # DR
       dr.form <- formula(paste(outcome.y[i.y],"~",
                                as.character(t.form[[2]]),"+",
                                as.character(x.form[2])))
       glm1 <- substitute(svyglm(formula=dr.form,design=sdesign.w,
                                 family=y.dist[i.y]))
       glm1 <- eval(glm1)
       results$glm.dr[[i.y]] <- glm1
       # transform from t to z for better FDR calculation
       results$z[i.y] <- 
          qnorm(pt(coef(summary(glm1))[2,"t value"],
                   glm1$df.residual))
   }
   if(verbose) cat("\nOutcome regression models complete\n")
### END OUTCOME ANALYSIS ###

### BEGIN TREATMENT CALCULATION ###
   results$effects <- vector("list",length(outcome.y))
   names(results$effects) <- attr(terms(y.form),"term.labels")

   # unweighted results
   sdesign   <- svydesign(ids=~1,weights=~samp.w,data=data0[i.treat,])
   means.un1 <- svymean(y.form,design=sdesign,na.rm=TRUE)
   sdesign   <- svydesign(ids=~1,weights=~samp.w,data=data0[!i.treat,])
   means.un0 <- svymean(y.form,design=sdesign,na.rm=TRUE)
   # propensity score results
   sdesign   <- svydesign(ids=~1,weights=~w,data=data0[!i.treat,])
   means.ps0 <- svymean(y.form,design=sdesign,na.rm=TRUE)

   for(i.y in 1:length(outcome.y))
   {
      results$effects[[i.y]] <-
         data.frame(E.y1 =rep(0,3),E.y0 =rep(0,3),
                    se.y1=rep(0,3),se.y0=rep(0,3),
                    TE   =rep(0,3),se.TE=rep(0,3))
      rownames(results$effects[[i.y]]) <- c("un","ps","dr")

      # collect unweighted statistics
      results$effects[[i.y]]$E.y1[]   <- means.un1[i.y]
      results$effects[[i.y]]$E.y0[1]  <- means.un0[i.y]
      results$effects[[i.y]]$se.y1[]  <- sqrt(diag(vcov(means.un1))[i.y])
      results$effects[[i.y]]$se.y0[1] <- sqrt(diag(vcov(means.un0))[i.y])

      results$effects[[i.y]]$TE[1]    <- coef(results$glm.un[[i.y]])[2]
      results$effects[[i.y]]$se.TE[1] <- sqrt(vcov(results$glm.un[[i.y]])[2,2])

      # collect PS statistics
      results$effects[[i.y]]$E.y0[2]  <- means.ps0[i.y]
      results$effects[[i.y]]$se.y0[2] <- sqrt(diag(vcov(means.ps0))[i.y])

      results$effects[[i.y]]$TE[2]    <- coef(results$glm.ps[[i.y]])[2]
      results$effects[[i.y]]$se.TE[2] <- sqrt(vcov(results$glm.ps[[i.y]])[2,2])

      # collect DR statistics
      results$effects[[i.y]]$TE[3]    <- coef(results$glm.dr[[i.y]])[2]
      results$effects[[i.y]]$se.TE[3] <- sqrt(vcov(results$glm.dr[[i.y]])[2,2])

      a <- data0[i.treat,]
      a[,as.character(t.form[2])] <- 0 # recode treated cases as control
      y.hat0 <- predict(results$glm.dr[[i.y]],
                        newdata=a,
                        type="response",
                        vcov=sum(i.treat)<1000) # use true Cov if Ntreat<1000
      results$effects[[i.y]]$E.y0[3]  <- mean(y.hat0)

      if(sum(i.treat)<1000) # use exact Cov for SE calculation
      {
         results$effects[[i.y]]$se.y0[3] <- sqrt(mean(vcov(y.hat0)))
      }
      else
      {
         # delta method? Tends to be too small
         # results$effects[[i.y]]$se.y0[3] <- sqrt(var(y.hat0)/length(y.hat0))
         # Monte Carlo the off-diagonals
         se.temp1 <- sum(vcov(y.hat0)) # vcov will return a vector here
         i <- sample(1:nrow(a), size=1000)
         y.hat0 <- predict(results$glm.dr[[i.y]],
                           newdata=a[i,],
                           type="response",
                           vcov=TRUE)
         se.temp2 <- sum(vcov(y.hat0))-sum(diag(vcov(y.hat0)))
         se.temp2 <- nrow(a)*(nrow(a)-1)*se.temp2/(length(i)*(length(i)-1))
         results$effects[[i.y]]$se.y0[3] <-
            sqrt((se.temp1+se.temp2)/(nrow(a)^2))
      }
   }

### END TREATMENT CALCULATION ###
   results$y.dist <- y.dist

   class(results) <- "fastDR"
   return(results)
}
