anova_1 <- function(exp_names,vcov_beta,beta, L3){
  L2 <- rep(0,length(exp_names))
  L2[L3] <- 1
  L <- L2
  VLbeta <- L %*% vcov_beta %*% L
  eig_VLbeta <- eigen(VLbeta)
  P <- eig_VLbeta$vectors
  d <- eig_VLbeta$values
  PtL <- crossprod(P, L)[1:1, ]
  CHI2 <- drop(PtL %*% beta)^2 / d[1:1]
  pvalue <- pchisq(q=CHI2, df=1, lower.tail=FALSE)
  data.frame("Df"=1, "Chisq"= CHI2, "Pr(>Chisq)"=pvalue)
  # }
}

#' Anova for a fitted \code{glmcat} model object
#' @description Compute an analysis of deviance table for one fitted \code{glmcat} model object.
#' @param object an object of class \code{"glmcat"}.
#' @param ... additional arguments.
#' @rdname anova
#' @method anova glmcat
#' @usage \method{anova}{glmcat}(object, ...)
#' @exportS3Method
anova.glmcat <-
  function(object, ...)
  {
    dots <- list(...)

    if(length(dots) == 0) {

      exp_names <- rownames(object$coef)
      beta <- coef(object, na.rm=TRUE)
      vcov_beta <- vcov(object)

      L1 <- 1:length(exp_names)
      L_in <- L1[-grep('Intercept',exp_names)]
      L_in <- as.list(L_in)

      table <- data.frame(do.call("rbind",(lapply(L_in, function(L) anova_1(exp_names,vcov_beta,beta,L)))))
      colnames(table) <- c("Df", "Chisq","Pr(>Chisq)")
      rownames(table) <- exp_names[unlist(L_in)]



      attr(table, "heading") <-
        "Wald test for glmcat models:\n"



    }else if(length(dots) == 1){

      mod2 <- dots[[1]]
      cat("Mod 1: ")
      print(object$formula)
      cat("Mod 2: ")
      print(mod2$formula)
      cat("\n")
      df = length(coef(mod2)) - length(coef(object))
      testStatistic <- 2 * (logLik(mod2) - logLik(object))
      mypval <- pchisq(testStatistic, df, lower.tail = FALSE)

      table <- data.frame("no.par"=c(object$`df of the model`,mod2$`df of the model`),
                 "AIC"= c(AIC(object),AIC(mod2)),
                 "logLik"= c(object$LogLikelihood,mod2$LogLikelihood),
                 "Chisq"=c(NA,testStatistic),
                 "Df" = c(NA,df),
                 "Pr(>Chisq)" = c(NA,mypval))

      colnames(table) <- c("no.par", "AIC","logLik","Chisq","Df","Pr(>Chisq)")
      rownames(table) <- c("Mod 1", "Mod 2")

      attr(table, "heading") <-
        "Likelihood ratio tests:\n"

    }
    class(table) <- c("anova.glmcat", "data.frame")

    return(table)
  }


#' Printing Anova for \code{glmcat} model fits
#' @description \code{print.anova} method for GLMcat objects.
#' @param x an object of class \code{"glmcat"}.
#' @param digits the number of digits in the printed table.
#' @param ... additional arguments affecting the summary produced.
#' @rdname print.anova
#' @method print anova.glmcat
#' @export
print.anova.glmcat <-
  function(x,
           digits=max(getOption("digits") - 2, 3),

           ...)
  {
    if (!is.null(heading <- attr(x, "heading")))
      cat(heading, "\n")

    signif.stars <- getOption("show.signif.stars")

    printCoefmat(x, digits=digits,
                 signif.stars=signif.stars,
                 # tst.ind=4, cs.ind=NULL, # zap.ind=2, #c(1,5),
                 P.values=TRUE, has.Pvalue=TRUE, na.print="", ...)
    return(invisible(x))
  }
